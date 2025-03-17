import React, { useState, useEffect, useRef, useMemo } from 'react';
import * as d3 from 'd3'; // d3 version 7.8+
import { Box, Typography, Skeleton, Slider, FormControl, FormLabel } from '@mui/material'; // @mui/material version 5.13+
import useWindowSize from '../../hooks/useWindowSize';
import { formatNumber, formatPropertyValue } from '../../utils/formatters';
import theme from '../../theme';
import Card from '../common/Card';
import { MoleculeProperty } from '../../types/molecule';

// Interface for Histogram bin data
interface HistogramBin {
  x0: number;
  x1: number;
  count: number;
  items: Array<any>;
}

// Props interface for the Histogram component
interface HistogramProps {
  data: Array<{ [key: string]: number | string }>;
  title: string;
  xAxisLabel: string;
  yAxisLabel: string;
  propertyName: string;
  propertyUnit?: string;
  color?: string;
  height?: number;
  width?: number | string;
  margin?: { top: number; right: number; bottom: number; left: number };
  showTooltip?: boolean;
  showGrid?: boolean;
  loading?: boolean;
  allowBinAdjustment?: boolean;
  initialBinCount?: number;
  minBinCount?: number;
  maxBinCount?: number;
  onBinClick?: (bin: any) => void;
  className?: string;
  style?: React.CSSProperties;
}

/**
 * A component that renders a histogram visualization for displaying the distribution
 * of a molecular property across a dataset of molecules.
 */
const Histogram: React.FC<HistogramProps> = ({
  data,
  title,
  xAxisLabel,
  yAxisLabel,
  propertyName,
  propertyUnit,
  color = theme.palette.primary.main,
  height = 300,
  width = '100%',
  margin = { top: 20, right: 20, bottom: 50, left: 50 },
  showTooltip = true,
  showGrid = true,
  loading = false,
  allowBinAdjustment = false,
  initialBinCount = 10,
  minBinCount = 5,
  maxBinCount = 30,
  onBinClick,
  className,
  style
}) => {
  // Refs for D3 manipulation
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  
  // State for bin count slider
  const [binCount, setBinCount] = useState<number>(initialBinCount);
  
  // Get window size for responsive adjustments
  const windowSize = useWindowSize();
  
  // Extract numerical values from data for histogram
  const processedData = useMemo(() => {
    if (!data || data.length === 0 || !propertyName) return [];
    
    return data
      .map(item => {
        const value = item[propertyName];
        // Ensure the value is a number
        return typeof value === 'number' ? value : 
               typeof value === 'string' ? parseFloat(value) : null;
      })
      .filter((value): value is number => value !== null && !isNaN(value));
  }, [data, propertyName]);
  
  // Generate histogram bins
  const histogramBins = useMemo(() => {
    if (processedData.length === 0) return [];
    
    try {
      // Calculate min and max values
      const minValue = Math.min(...processedData);
      const maxValue = Math.max(...processedData);
      
      // Handle edge case where all values are the same
      if (minValue === maxValue) {
        return [{
          x0: minValue - 0.5,
          x1: maxValue + 0.5,
          count: processedData.length,
          items: [...data]
        }];
      }
      
      // Create histogram generator with specified bin count
      const binGenerator = d3.bin<number, number>()
        .domain([minValue, maxValue])
        .thresholds(Array.from({ length: binCount - 1 }, (_, i) => 
          minValue + ((maxValue - minValue) / binCount) * (i + 1)
        ));
      
      const bins = binGenerator(processedData);
      
      // Enhance bins with original items
      return bins.map(bin => {
        const binStart = bin.x0 ?? 0;
        const binEnd = bin.x1 ?? 0;
        
        const items = data.filter(d => {
          const value = d[propertyName];
          const numValue = typeof value === 'number' ? value : 
                         typeof value === 'string' ? parseFloat(value) : null;
          return numValue !== null && !isNaN(numValue) && 
                 numValue >= binStart && numValue < binEnd;
        });
        
        return {
          x0: binStart,
          x1: binEnd,
          count: bin.length,
          items
        };
      });
    } catch (error) {
      console.error('Error generating histogram bins:', error);
      return [];
    }
  }, [processedData, binCount, data, propertyName]);
  
  // Handle bin count slider change
  const handleBinCountChange = (_event: React.ChangeEvent<HTMLInputElement>, value: number | number[]) => {
    if (typeof value === 'number') {
      setBinCount(value);
    }
  };
  
  // Create and update the histogram visualization
  useEffect(() => {
    if (!svgRef.current || histogramBins.length === 0) return;
    
    try {
      // Clear existing SVG content
      const svg = d3.select(svgRef.current);
      svg.selectAll('*').remove();
      
      // Calculate dimensions
      let chartWidth: number;
      if (typeof width === 'string') {
        if (width.includes('%')) {
          const percentValue = parseInt(width);
          chartWidth = svgRef.current.parentElement 
            ? (svgRef.current.parentElement.clientWidth * percentValue / 100)
            : 500;
        } else {
          chartWidth = parseFloat(width);
        }
      } else {
        chartWidth = width || 500;
      }
      
      const chartHeight = height;
      
      // Set SVG dimensions
      svg.attr('width', chartWidth)
         .attr('height', chartHeight);
      
      // Create chart group with margins
      const chart = svg.append('g')
        .attr('transform', `translate(${margin.left},${margin.top})`);
      
      // Calculate inner dimensions
      const innerWidth = chartWidth - margin.left - margin.right;
      const innerHeight = chartHeight - margin.top - margin.bottom;
      
      // Find data extent for x and y scales
      const xMin = d3.min(histogramBins, d => d.x0) || 0;
      const xMax = d3.max(histogramBins, d => d.x1) || 1;
      const maxCount = d3.max(histogramBins, d => d.count) || 0;
      
      // Create x scale
      const xScale = d3.scaleLinear()
        .domain([xMin, xMax])
        .range([0, innerWidth])
        .nice();
      
      // Create y scale
      const yScale = d3.scaleLinear()
        .domain([0, maxCount * 1.1]) // Add 10% padding at top
        .range([innerHeight, 0])
        .nice();
      
      // Add x axis
      chart.append('g')
        .attr('class', 'x-axis')
        .attr('transform', `translate(0,${innerHeight})`)
        .call(
          d3.axisBottom(xScale)
            .ticks(5)
            .tickFormat(d => formatNumber(+d))
        );
      
      // Add x axis label
      chart.append('text')
        .attr('class', 'x-axis-label')
        .attr('x', innerWidth / 2)
        .attr('y', innerHeight + 40)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .text(xAxisLabel);
      
      // Add y axis
      chart.append('g')
        .attr('class', 'y-axis')
        .call(
          d3.axisLeft(yScale)
            .ticks(5)
            .tickFormat(d => formatNumber(+d))
        );
      
      // Add y axis label
      chart.append('text')
        .attr('class', 'y-axis-label')
        .attr('transform', 'rotate(-90)')
        .attr('x', -innerHeight / 2)
        .attr('y', -40)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .text(yAxisLabel);
      
      // Add grid lines if enabled
      if (showGrid) {
        // Add x grid lines
        chart.append('g')
          .attr('class', 'grid x-grid')
          .attr('transform', `translate(0,${innerHeight})`)
          .call(
            d3.axisBottom(xScale)
              .tickSize(-innerHeight)
              .tickFormat(() => '')
          )
          .selectAll('line')
          .style('stroke', theme.palette.divider)
          .style('stroke-opacity', 0.5);
        
        // Add y grid lines
        chart.append('g')
          .attr('class', 'grid y-grid')
          .call(
            d3.axisLeft(yScale)
              .tickSize(-innerWidth)
              .tickFormat(() => '')
          )
          .selectAll('line')
          .style('stroke', theme.palette.divider)
          .style('stroke-opacity', 0.5);
      }
      
      // Create and render bars
      const bars = chart.selectAll('.bar')
        .data(histogramBins)
        .enter()
        .append('rect')
        .attr('class', 'bar')
        .attr('x', d => xScale(d.x0))
        .attr('y', d => yScale(0)) // Start at bottom for animation
        .attr('width', d => Math.max(0, xScale(d.x1) - xScale(d.x0) - 1))
        .attr('height', 0) // Start with height 0 for animation
        .attr('fill', color)
        .style('cursor', onBinClick ? 'pointer' : 'default');
      
      // Add event listeners
      if (showTooltip) {
        bars.on('mouseover', function(event: MouseEvent, d: HistogramBin) {
          // Change bar color on hover
          d3.select(this).attr('fill', d3.color(color)!.darker(0.2).toString());
          
          // Show tooltip
          if (tooltipRef.current) {
            const tooltip = tooltipRef.current;
            
            // Format bin range for display
            const minValue = formatPropertyValue(d.x0, propertyName, propertyUnit);
            const maxValue = formatPropertyValue(d.x1, propertyName, propertyUnit);
            
            tooltip.innerHTML = `
              <div>Range: ${minValue} - ${maxValue}</div>
              <div>Count: ${d.count} molecule${d.count !== 1 ? 's' : ''}</div>
            `;
            
            tooltip.style.display = 'block';
            tooltip.style.left = `${event.pageX + 10}px`;
            tooltip.style.top = `${event.pageY - 40}px`;
          }
        })
        .on('mouseout', function() {
          // Restore original bar color
          d3.select(this).attr('fill', color);
          
          // Hide tooltip
          if (tooltipRef.current) {
            tooltipRef.current.style.display = 'none';
          }
        });
      }
      
      // Add click handler
      if (onBinClick) {
        bars.on('click', function(_event: MouseEvent, d: HistogramBin) {
          onBinClick(d);
        });
      }
      
      // Animate bars
      bars.transition()
        .duration(500)
        .attr('y', d => yScale(d.count))
        .attr('height', d => innerHeight - yScale(d.count));
    } catch (error) {
      console.error('Error rendering histogram:', error);
    }
    
    // Cleanup function to remove event listeners
    return () => {
      if (svgRef.current) {
        d3.select(svgRef.current).selectAll('*').remove();
      }
    };
  }, [histogramBins, data, propertyName, propertyUnit, width, height, margin, showGrid, xAxisLabel, yAxisLabel, color, onBinClick, showTooltip, windowSize]);
  
  return (
    <Card className={className} style={style}>
      <Typography variant="h6" gutterBottom>
        {title}
      </Typography>
      
      {allowBinAdjustment && (
        <Box sx={{ mb: 2 }}>
          <FormControl fullWidth size="small">
            <FormLabel>Bin Count: {binCount}</FormLabel>
            <Slider
              value={binCount}
              min={minBinCount}
              max={maxBinCount}
              step={1}
              onChange={handleBinCountChange}
              aria-labelledby="bin-count-slider"
              size="small"
            />
          </FormControl>
        </Box>
      )}
      
      {loading ? (
        <Skeleton variant="rectangular" width="100%" height={height} />
      ) : processedData.length === 0 ? (
        <Box 
          sx={{ 
            height, 
            width: '100%', 
            display: 'flex', 
            alignItems: 'center', 
            justifyContent: 'center',
            color: 'text.secondary'
          }}
        >
          <Typography variant="body2">
            No numerical data available for this property
          </Typography>
        </Box>
      ) : (
        <Box sx={{ position: 'relative', width: '100%', height: height }}>
          <svg ref={svgRef} width="100%" height={height} />
          <div 
            ref={tooltipRef} 
            style={{
              display: 'none',
              position: 'fixed',
              backgroundColor: theme.palette.background.paper,
              padding: '8px',
              borderRadius: '4px',
              boxShadow: theme.shadows[2],
              pointerEvents: 'none',
              zIndex: 1000,
              fontSize: '12px'
            }}
          />
        </Box>
      )}
    </Card>
  );
};

export default Histogram;