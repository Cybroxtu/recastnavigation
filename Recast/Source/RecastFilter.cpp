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
#include "RecastAssert.h"

/// @par
///
/// Allows the formation of walkable regions that will flow over low lying 
/// objects such as curbs, and up structures such as stairways. 
/// 
/// Two neighboring spans are walkable if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) < waklableClimb</tt>
/// 
/// @warning Will override the effect of #rcFilterLedgeSpans.  So if both filters are used, call
/// #rcFilterLedgeSpans after calling this filter. 
///
/// @see rcHeightfield, rcConfig
void rcFilterLowHangingWalkableObstacles(rcContext* ctx, const int walkableClimb, rcHeightfield& solid)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_FILTER_LOW_OBSTACLES);
	
	const int w = solid.width;
	const int h = solid.height;
	
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			rcSpan* ps = 0;
			bool previousWalkable = false;
			unsigned char previousArea = RC_NULL_AREA;
			
			for (rcSpan* s = solid.spans[x + y*w]; s; ps = s, s = s->next)
			{
				const bool walkable = s->area != RC_NULL_AREA;
				// If current span is not walkable, but there is walkable
				// span just below it, mark the span above it walkable too.
				if (!walkable && previousWalkable)
				{
					if (rcAbs((int)s->smax - (int)ps->smax) <= walkableClimb)
						s->area = previousArea;
				}
				// Copy walkable flag so that it cannot propagate
				// past multiple non-walkable objects.
				previousWalkable = walkable;
				previousArea = s->area;
			}
		}
	}
}

/// @par
///
/// A ledge is a span with one or more neighbors whose maximum is further away than @p walkableClimb
/// from the current span's maximum.
/// This method removes the impact of the overestimation of conservative voxelization 
/// so the resulting mesh will not have regions hanging in the air over ledges.
/// 
/// A span is a ledge if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) > walkableClimb</tt>
/// 
/// @see rcHeightfield, rcConfig
void rcFilterLedgeSpans(rcContext* ctx, const int walkableHeight, const int walkableClimb,
						rcHeightfield& solid)
{
	rcAssert(ctx);
	
	rcScopedTimer timer(ctx, RC_TIMER_FILTER_BORDER);

	const int w = solid.width;
	const int h = solid.height;
	const int MAX_HEIGHT = 0xffff;
	
	// Mark border spans.
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
		    // 遍历所有的 solid span
			for (rcSpan* s = solid.spans[x + y*w]; s; s = s->next)
			{
				// Skip non walkable spans.
				if (s->area == RC_NULL_AREA)
					continue;

				// 这里是 open span 的底部、顶部高度
				const int bot = (int)(s->smax);
				const int top = s->next ? (int)(s->next->smin) : MAX_HEIGHT;
				
				// Find neighbours minimum height.
				// 邻接 open span 的底部高度 - 当前 open span 的底部高度
				int minh = MAX_HEIGHT;

				// Min and max height of accessible neighbours.
				int asmin = s->smax;
				int asmax = s->smax;

				// 遍历邻接四方向
				for (int dir = 0; dir < 4; ++dir)
				{
					int dx = x + rcGetDirOffsetX(dir);
					int dy = y + rcGetDirOffsetY(dir);
					// Skip neighbours which are out of bounds.
					if (dx < 0 || dy < 0 || dx >= w || dy >= h)
					{
						minh = rcMin(minh, -walkableClimb - bot);
						continue;
					}

					// From minus infinity to the first span.
					rcSpan* ns = solid.spans[dx + dy*w];
                    // 邻接 open span 的底部、顶部高度，第一个 open span 的高度从 负无穷大到第一个 solid span 的底部高度
					int nbot = -walkableClimb; // 说是负无穷大，其实小于 -walkableClimb 的值没有意义，直接使用 -walkableClimb 了
					int ntop = ns ? (int)ns->smin : MAX_HEIGHT;
					// Skip neightbour if the gap between the spans is too small.
					// 判断两个 open span 在 y 轴上相交的部分高度，如果小于 walkableHeight，说明这两个 open span 间不能通过
					if (rcMin(top,ntop) - rcMax(bot,nbot) > walkableHeight)
						minh = rcMin(minh, nbot - bot);
					
					// Rest of the spans.
					for (ns = solid.spans[dx + dy*w]; ns; ns = ns->next)
					{
						nbot = (int)ns->smax;
						ntop = ns->next ? (int)ns->next->smin : MAX_HEIGHT;
						// Skip neightbour if the gap between the spans is too small.
						if (rcMin(top,ntop) - rcMax(bot,nbot) > walkableHeight)
						{
							minh = rcMin(minh, nbot - bot);
						
							// Find min/max accessible neighbour height. 
							if (rcAbs(nbot - bot) <= walkableClimb)
							{
							    // 如果邻接 open span 在当前 open span 的高度差 walkableClimb 以内
							    // 则更新可访问的邻接 span 最小、最大坐标值
								if (nbot < asmin) asmin = nbot;
								if (nbot > asmax) asmax = nbot;
							}
							
						}
					}
				}
				
				// The current span is close to a ledge if the drop to any
				// neighbour span is less than the walkableClimb.
				if (minh < -walkableClimb)
				{
                    // 因为是邻接高度减去当前高度，所以这个判断代表邻接 span 比当前 span 要低
                    // 此时将当前 span 标记为不可行走，也就是说一高一低两个 span，按定义两个都是 ledge
                    // 但是只将较高的 span 标记为不可行走
                    // 这里主要是为了把悬崖边缘给标记为不可行走区域
					s->area = RC_NULL_AREA;
				}
				// If the difference between all neighbours is too large,
				// we are at steep slope, mark the span as ledge.
				else if ((asmax - asmin) > walkableClimb)
				{
				    // 如果可达的邻接的最大、最小高度的插值大于了 walkableClimb
				    // 则代表当前是在一个陡坡的中间，应该将当前 span 标记为不可行走
				    // 这里是为了解决之前 rcFilterLowHangingWalkableObstacles 带来的问题
				    // 如果一整个陡坡被连带着标记为了可以行走，这里需要进行纠正
					s->area = RC_NULL_AREA;
				}
			}
		}
	}
}

/// @par
///
/// For this filter, the clearance above the span is the distance from the span's 
/// maximum to the next higher span's minimum. (Same grid column.)
/// 
/// @see rcHeightfield, rcConfig
void rcFilterWalkableLowHeightSpans(rcContext* ctx, int walkableHeight, rcHeightfield& solid)
{
	rcAssert(ctx);
	
	rcScopedTimer timer(ctx, RC_TIMER_FILTER_WALKABLE);
	
	const int w = solid.width;
	const int h = solid.height;
	const int MAX_HEIGHT = 0xffff;
	
	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			for (rcSpan* s = solid.spans[x + y*w]; s; s = s->next)
			{
				const int bot = (int)(s->smax);
				const int top = s->next ? (int)(s->next->smin) : MAX_HEIGHT;
				if ((top - bot) <= walkableHeight)
					s->area = RC_NULL_AREA;
			}
		}
	}
}
