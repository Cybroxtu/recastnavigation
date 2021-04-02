//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"

// 判断两个 AABB 包围盒是否相交
inline bool overlapBounds(const float* amin, const float* amax, const float* bmin, const float* bmax)
{
	bool overlap = true;
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
	return overlap;
}

// 判断一维坐标轴上两条线段是否相交
inline bool overlapInterval(unsigned short amin, unsigned short amax,
							unsigned short bmin, unsigned short bmax)
{
	if (amax < bmin) return false;
	if (amin > bmax) return false;
	return true;
}


static rcSpan* allocSpan(rcHeightfield& hf)
{
	// If running out of memory, allocate new page and update the freelist.
	if (!hf.freelist || !hf.freelist->next)
	{
		// Create new page.
		// Allocate memory for the new pool.
		rcSpanPool* pool = (rcSpanPool*)rcAlloc(sizeof(rcSpanPool), RC_ALLOC_PERM);
		if (!pool) return 0;

		// Add the pool into the list of pools.
		pool->next = hf.pools;
		hf.pools = pool;
		// Add new items to the free list.
		rcSpan* freelist = hf.freelist;
		rcSpan* head = &pool->items[0];
		rcSpan* it = &pool->items[RC_SPANS_PER_POOL];
		do
		{
			--it;
			it->next = freelist;
			freelist = it;
		}
		while (it != head);
		hf.freelist = it;
	}
	
	// Pop item from in front of the free list.
	rcSpan* it = hf.freelist;
	hf.freelist = hf.freelist->next;
	return it;
}

static void freeSpan(rcHeightfield& hf, rcSpan* ptr)
{
	if (!ptr) return;
	// Add the node in front of the free list.
	ptr->next = hf.freelist;
	hf.freelist = ptr;
}

// 将新的 solid span 添加到高度场内
static bool addSpan(rcHeightfield& hf, const int x, const int y,
					const unsigned short smin, const unsigned short smax,
					const unsigned char area, const int flagMergeThr)
{
	int idx = x + y*hf.width;
	
	rcSpan* s = allocSpan(hf);
	if (!s)
		return false;
	s->smin = smin;
	s->smax = smax;
	s->area = area;
	s->next = 0;
	
	// Empty cell, add the first span.
	if (!hf.spans[idx])
	{
		hf.spans[idx] = s;
		return true;
	}
	rcSpan* prev = 0;
	rcSpan* cur = hf.spans[idx];
	
	// Insert and merge spans.
	// 将新添加的 span 根据其高度插入到链表中
    // 需要注意的事，如果新添加的 span 与老的 span 存在重叠，将老的 span 合并到新 span 里
    // 1. 如果两个 span 的顶部高度差在合并范围 flagMergeThr 以内，则使用两个 span 中 area 较大的值（其实就是 RC_WALKABLE_AREA 了）
    // 2. 如果两个 span 的顶部高度差在合并范围 flagMergeThr 以外，则使用新增 span 的 area 值
	while (cur)
	{
		if (cur->smin > s->smax)
		{
			// Current span is further than the new span, break.
			break;
		}
		else if (cur->smax < s->smin)
		{
			// Current span is before the new span advance.
			prev = cur;
			cur = cur->next;
		}
		else
		{
			// Merge spans.
			if (cur->smin < s->smin)
				s->smin = cur->smin;
			if (cur->smax > s->smax)
				s->smax = cur->smax;
			
			// Merge flags.
			if (rcAbs((int)s->smax - (int)cur->smax) <= flagMergeThr)
				s->area = rcMax(s->area, cur->area);
			
			// Remove current span.
			rcSpan* next = cur->next;
			freeSpan(hf, cur);
			if (prev)
				prev->next = next;
			else
				hf.spans[idx] = next;
			cur = next;
		}
	}
	
	// Insert new span.
	if (prev)
	{
		s->next = prev->next;
		prev->next = s;
	}
	else
	{
		s->next = hf.spans[idx];
		hf.spans[idx] = s;
	}

	return true;
}

/// @par
///
/// The span addition can be set to favor flags. If the span is merged to
/// another span and the new @p smax is within @p flagMergeThr units
/// from the existing span, the span flags are merged.
///
/// @see rcHeightfield, rcSpan.
bool rcAddSpan(rcContext* ctx, rcHeightfield& hf, const int x, const int y,
			   const unsigned short smin, const unsigned short smax,
			   const unsigned char area, const int flagMergeThr)
{
	rcAssert(ctx);

	if (!addSpan(hf, x, y, smin, smax, area, flagMergeThr))
	{
		ctx->log(RC_LOG_ERROR, "rcAddSpan: Out of memory.");
		return false;
	}

	return true;
}

// divides a convex polygons into two convex polygons on both sides of a line
// 将一个凸多边形沿着一条轴线切割成两个凸多边形
// 参数：
// in - 输入凸多边形的顶点数组
// nin - 输入凸多边形的顶点数量
// out1 - 输出凸多边形 1 的顶点数组
// nout1 - 输出凸多边形 1 的顶点数量
// out2 - 输出凸多边形 2 的顶点数组
// nout2 - 输出凸多边形 2 的顶点数量
// x - 要进行切割的轴坐标
// axis - 代表要在哪个坐标轴上进行切割，顶点坐标数组的偏移量，0-x 1-y 2-z
static void dividePoly(const float* in, int nin,
					  float* out1, int* nout1,
					  float* out2, int* nout2,
					  float x, int axis)
{
    // 计算出多边形各顶点坐标在 axis 轴上与 x 的差值
    // 用于判断两个点在 axis 轴上是否相等，是的话，说明这两个点一定在切割后的同一边
	float d[12];
	for (int i = 0; i < nin; ++i)
		d[i] = x - in[i*3+axis];

	// m - 切割出的多边形的顶点数量
    // n - 剩下的多边形的顶点数量
	int m = 0, n = 0;

    // 遍历多边形的每一条边 j-i
    // 每次遍历一条边时，除了新添加的割点外，只处理 j-i 的结束点 i
    for (int i = 0, j = nin-1; i < nin; j=i, ++i)
	{
		bool ina = d[j] >= 0;
		bool inb = d[i] >= 0;
        // 这里判断两个点是否在切割轴线 x 的同一边
		if (ina != inb)
		{
            // 两个顶点在切割线的不同边，代表可以切割
            // 求出 j-i 与切割线的交点，这个交点被切割后的两个多边形所共享
			float s = d[j] / (d[j] - d[i]);
			out1[m*3+0] = in[j*3+0] + (in[i*3+0] - in[j*3+0])*s;
			out1[m*3+1] = in[j*3+1] + (in[i*3+1] - in[j*3+1])*s;
			out1[m*3+2] = in[j*3+2] + (in[i*3+2] - in[j*3+2])*s;
			rcVcopy(out2 + n*3, out1 + m*3); // 这个
			m++;
			n++;

			// add the i'th point to the right polygon. Do NOT add points that are on the dividing line
			// since these were already added above
			if (d[i] > 0)
			{
                // 点 i 在分割线左边
				rcVcopy(out1 + m*3, in + i*3);
				m++;
			}
			else if (d[i] < 0)
			{
                // 点 i 在分割线右边
				rcVcopy(out2 + n*3, in + i*3);
				n++;
			}
		}
		else // same side
		{
		    // 两个顶点在切割线的同一边，包括在分割线上

			// add the i'th point to the right polygon. Addition is done even for points on the dividing line
			if (d[i] >= 0)
			{
                // 1. 点 i 在分割线左边，此时它只属于第一个多边形，然后可以进行下一条边的遍历
                // 2. 点 i 在分割线上，此时它同时属于两个多边形
				rcVcopy(out1 + m*3, in + i*3);
				m++;
				if (d[i] != 0) // 两个顶点都在分割线上
					continue;
			}

			// 1. 两个点都在分割线右边
			// 此时第一个多边形没有可添加的点，第二个多边形将每条边 j-i 的结束点添加到多边形内
			// 2. 两个点都在分割线上，此将边 j-i 的结束点添加到第二个多边形内
			rcVcopy(out2 + n*3, in + i*3);
			n++;
		}
	}

	*nout1 = m;
	*nout2 = n;
}



static bool rasterizeTri(const float* v0, const float* v1, const float* v2,
						 const unsigned char area, rcHeightfield& hf,
						 const float* bmin, const float* bmax,
						 const float cs, const float ics, const float ich,
						 const int flagMergeThr)
{
	const int w = hf.width;
	const int h = hf.height;
	float tmin[3], tmax[3];
	const float by = bmax[1] - bmin[1];
	
	// Calculate the bounding box of the triangle.
	// 计算出三角形的 AABB 包围盒
	rcVcopy(tmin, v0);
	rcVcopy(tmax, v0);
	rcVmin(tmin, v1);
	rcVmin(tmin, v2);
	rcVmax(tmax, v1);
	rcVmax(tmax, v2);
	
	// If the triangle does not touch the bbox of the heightfield, skip the triagle.
	// 判断三角形的包围盒是否与高度场的包围盒相交
	// 如果不相交，则表示这个三角形不在当前构建导航网格的空间内，可以直接跳过
	if (!overlapBounds(bmin, bmax, tmin, tmax))
		return true;
	
	// Calculate the footprint of the triangle on the grid's y-axis
	// ics 是 cellSize 的倒数，用于将原始坐标转换为高度场的格子坐标
	// 这里 y0 y1 分别是三角形包围盒与高度场的包围盒最小点在 z 轴坐标上的差值转换到格子坐标后的值
	// 方便后面计算
	int y0 = (int)((tmin[2] - bmin[2])*ics);
	int y1 = (int)((tmax[2] - bmin[2])*ics);
	y0 = rcClamp(y0, 0, h-1);
	y1 = rcClamp(y1, 0, h-1);
	
	// Clip the triangle into all grid cells it touches.
	float buf[7*3*4];
	float *in = buf, *inrow = buf+7*3, *p1 = inrow+7*3, *p2 = p1+7*3;

	// 准备初始的三个顶点数据
	rcVcopy(&in[0], v0);
	rcVcopy(&in[1*3], v1);
	rcVcopy(&in[2*3], v2);
	int nvrow, nvIn = 3;

	// 遍历三角形包围盒内的格子，在 z 轴上进行切割
	for (int y = y0; y <= y1; ++y)
	{
		// Clip polygon to row. Store the remaining polygon as well
		const float cz = bmin[2] + y*cs; // 三角形包围盒的起始 z 轴值 + 当前格子偏移值

		// 对 in/nvIn 代表的三角形进行切割
		// cz+cs 代表每次切割下一个体素大小的区域
		// axis = 2 代表在 z 轴上进行切割
		// 切割出来要进行光栅化的多边形为 inrow/nvrow
        // 待继续切割的多边形为 p1/nvIn
        // 所以这里将 in 与 p1 进行交换
		dividePoly(in, nvIn, inrow, &nvrow, p1, &nvIn, cz+cs, 2);
		rcSwap(in, p1);
		if (nvrow < 3) continue; // 没有割到东西

		// find the horizontal bounds in the row
		float minX = inrow[0], maxX = inrow[0];
		for (int i=1; i<nvrow; ++i)
		{
			if (minX > inrow[i*3])	minX = inrow[i*3];
			if (maxX < inrow[i*3])	maxX = inrow[i*3];
		}

		// 计算出切割出的待光栅化多边形里，其包围盒与高度场包围盒在 x 轴上的差值
		int x0 = (int)((minX - bmin[0])*ics);
		int x1 = (int)((maxX - bmin[0])*ics);
		x0 = rcClamp(x0, 0, w-1);
		x1 = rcClamp(x1, 0, w-1);

		int nv, nv2 = nvrow;

        // 遍历切割出的多边形，在 x 轴上对其再次进行切割
		for (int x = x0; x <= x1; ++x)
		{
			// Clip polygon to column. store the remaining polygon as well
			const float cx = bmin[0] + x*cs;
			dividePoly(inrow, nv2, p1, &nv, p2, &nv2, cx+cs, 0);
			rcSwap(inrow, p2);
			if (nv < 3) continue;
			
			// Calculate min and max of the span.
			float smin = p1[1], smax = p1[1];
			for (int i = 1; i < nv; ++i)
			{
				smin = rcMin(smin, p1[i*3+1]);
				smax = rcMax(smax, p1[i*3+1]);
			}
			smin -= bmin[1];
			smax -= bmin[1];
			// Skip the span if it is outside the heightfield bbox
			if (smax < 0.0f) continue;
			if (smin > by) continue;
			// Clamp the span to the heightfield bbox.
			if (smin < 0.0f) smin = 0;
			if (smax > by) smax = by;
			
			// Snap the span to the heightfield height grid.
			unsigned short ismin = (unsigned short)rcClamp((int)floorf(smin * ich), 0, RC_SPAN_MAX_HEIGHT);
			unsigned short ismax = (unsigned short)rcClamp((int)ceilf(smax * ich), (int)ismin+1, RC_SPAN_MAX_HEIGHT);

			// 到这里，代表三角形面与对应格子相交，应该将相关格子的数据添加到高度场中
			if (!addSpan(hf, x, y, ismin, ismax, area, flagMergeThr))
				return false;
		}
	}

	return true;
}

/// @par
///
/// No spans will be added if the triangle does not overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangle(rcContext* ctx, const float* v0, const float* v1, const float* v2,
						 const unsigned char area, rcHeightfield& solid,
						 const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);

	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	if (!rasterizeTri(v0, v1, v2, area, solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
	{
		ctx->log(RC_LOG_ERROR, "rcRasterizeTriangle: Out of memory.");
		return false;
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const int /*nv*/,
						  const int* tris, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	// Rasterize triangles.
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[tris[i*3+0]*3];
		const float* v1 = &verts[tris[i*3+1]*3];
		const float* v2 = &verts[tris[i*3+2]*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const int /*nv*/,
						  const unsigned short* tris, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;

	// Rasterize triangles.
	// 遍历所有三角形，进行光栅化处理
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[tris[i*3+0]*3];
		const float* v1 = &verts[tris[i*3+1]*3];
		const float* v2 = &verts[tris[i*3+2]*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);
	
	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	// Rasterize triangles.
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[(i*3+0)*3];
		const float* v1 = &verts[(i*3+1)*3];
		const float* v2 = &verts[(i*3+2)*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}
