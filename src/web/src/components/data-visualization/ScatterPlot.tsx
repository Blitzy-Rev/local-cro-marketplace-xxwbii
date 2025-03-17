import React, { useState, useEffect, useRef, useMemo } from 'react';
import * as d3 from 'd3'; // version 7.8+
import { Box, Typography, Tooltip, Skeleton, FormControl, FormLabel, Switch, FormControlLabel } from '@mui/material'; // v5.13+
import useWindowSize from '../../hooks/useWindowSize';
import { formatNumber, formatPropertyValue } from '../../utils/formatters';
import { getPropertyDisplayName, getPropertyUnit } from '../../utils/molecularUtils';
import theme from '../../theme';
import Card from '../common/Card';

/**
 * Props interface for the ScatterPlot component
 */
interface ScatterPlotProps {
  /** Array of data objects containing property values to plot */
  data: Array<{ [key: string]: number | string }>;
  /** Title of the scatter plot */
  title?: string;
  /** Property name to use for x-axis */
  xAxisProperty: string;
  /** Property name to use for y-axis */
  yAxisProperty: string;
  /** Label for x-axis (defaults to formatted property name) */
  xAxisLabel?: string;
  /** Label for y-axis (defaults to formatted property name) */
  yAxisLabel?: string;
  /** Color for data points (defaults to primary color) */
  color?: string;
  /** Height of the chart in pixels */
  height?: number;
  /** Width of the chart (can be number or percentage string) */
  width?: number | string;
  /** Margins around the chart */
  margin?: { top: number; right: number; bottom: number; left: number };
  /** Whether to show tooltips on hover */
  showTooltip?: boolean;
  /** Whether to show grid lines */
  showGrid?: boolean;
  /** Whether the data is still loading */
  loading?: boolean;
  /** Whether to allow zooming and panning */
  allowZoom?: boolean;
  /** Callback when a point is clicked */
  onPointClick?: (point: any) => void;
  /** Additional class name */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
}

/**
 * Interface representing a data point in the scatter plot
 */
interface ScatterPoint {
  /** X coordinate (derived from xAxisProperty) */
  x: number;
  /** Y coordinate (derived from yAxisProperty) */
  y: number;
  /** Label for the point (usually SMILES or name) */
  label: string;
  /** Reference to the original data object */
  originalData: any;
}

/**
 * A component that renders a scatter plot visualization for displaying 
 * relationships between two molecular properties.
 * 
 * This component provides interactive features like tooltips, point selection,
 * zooming, and customizable styling to help researchers analyze correlations
 * between different molecular properties.
 */
const ScatterPlot: React.FC<ScatterPlotProps> = ({
  data,
  title,
  xAxisProperty,
  yAxisProperty,
  xAxisLabel,
  yAxisLabel,
  color = theme.palette.primary.main,
  height = 300,
  width = '100%',
  margin = { top: 20, right: 20, bottom: 50, left: 50 },
  showTooltip = true,
  showGrid = true,
  loading = false,
  allowZoom = false,
  onPointClick,
  className,
  style,
}) => {
  // Refs for accessing DOM elements
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  
  // Get window dimensions for responsive sizing
  const windowSize = useWindowSize();
  
  // State for zoom transform
  const [zoomState, setZoomState] = useState({
    scale: 1,
    translateX: 0,
    translateY: 0
  });
  
  // Process data into a format suitable for D3
  const processedData = useMemo<ScatterPoint[]>(() => {
    if (!data || data.length === 0) return [];
    
    return data
      .filter(item => {
        // Ensure both properties exist and are numeric
        const xValue = Number(item[xAxisProperty]);
        const yValue = Number(item[yAxisProperty]);
        return !isNaN(xValue) && !isNaN(yValue);
      })
      .map(item => {
        // Create a ScatterPoint object for each data item
        const xValue = Number(item[xAxisProperty]);
        const yValue = Number(item[yAxisProperty]);
        
        // Use SMILES or a property named 'name' or 'id' for the label
        const label = (item.smiles || item.name || item.id || '').toString();
        
        return {
          x: xValue,
          y: yValue,
          label,
          originalData: item
        };
      });
  }, [data, xAxisProperty, yAxisProperty]);
  
  // Handle click events on data points
  const handlePointClick = (point: ScatterPoint, event: React.MouseEvent) => {
    if (onPointClick) {
      onPointClick(point.originalData);
    }
  };
  
  // Display tooltip with data point information
  const showTooltipFunc = (point: ScatterPoint, event: MouseEvent) => {
    if (!showTooltip || !tooltipRef.current) return;
    
    const tooltip = tooltipRef.current;
    
    // Position tooltip near mouse pointer
    const x = event.clientX;
    const y = event.clientY;
    
    tooltip.style.left = `${x + 10}px`;
    tooltip.style.top = `${y - 10}px`;
    tooltip.style.display = 'block';
    
    // Format x and y values with appropriate units
    const xUnit = getPropertyUnit(xAxisProperty);
    const yUnit = getPropertyUnit(yAxisProperty);
    
    const xDisplayName = getPropertyDisplayName(xAxisProperty);
    const yDisplayName = getPropertyDisplayName(yAxisProperty);
    
    const xValueFormatted = formatPropertyValue(point.x, xAxisProperty, xUnit);
    const yValueFormatted = formatPropertyValue(point.y, yAxisProperty, yUnit);
    
    // Set tooltip content
    tooltip.innerHTML = `
      <div style="padding: 8px;">
        <div style="font-weight: bold;">${point.label || 'Molecule'}</div>
        <div>${xDisplayName}: ${xValueFormatted}</div>
        <div>${yDisplayName}: ${yValueFormatted}</div>
      </div>
    `;
  };
  
  // Hide the tooltip
  const hideTooltip = () => {
    if (tooltipRef.current) {
      tooltipRef.current.style.display = 'none';
    }
  };
  
  // Handle zoom and pan events on the scatter plot
  const handleZoom = (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
    if (!allowZoom) return;
    
    const transform = event.transform;
    
    setZoomState({
      scale: transform.k,
      translateX: transform.x,
      translateY: transform.y
    });
    
    // Apply transform to the chart group
    d3.select(svgRef.current)
      .select('g.chart-group')
      .attr('transform', transform.toString());
  };
  
  // Reset zoom and pan to default state
  const resetZoom = () => {
    if (!allowZoom || !svgRef.current) return;
    
    setZoomState({
      scale: 1,
      translateX: 0,
      translateY: 0
    });
    
    d3.select(svgRef.current)
      .transition()
      .duration(750)
      .call(
        // @ts-ignore (zoom.transform expects d3.ZoomBehavior which is hard to type correctly)
        d3.zoom().transform,
        d3.zoomIdentity
      );
  };
  
  // Create and update the scatter plot using D3
  useEffect(() => {
    if (!svgRef.current) return;
    
    const svg = d3.select(svgRef.current);
    
    // Clear previous content
    svg.selectAll('*').remove();
    
    // If no data or processing error, show empty message
    if (!processedData || processedData.length === 0) {
      svg.append('text')
        .attr('x', '50%')
        .attr('y', '50%')
        .attr('text-anchor', 'middle')
        .attr('dominant-baseline', 'middle')
        .style('font-size', '14px')
        .style('fill', theme.palette.text.secondary)
        .text('No data available');
      return;
    }
    
    // Calculate the actual dimensions of the chart
    const containerWidth = typeof width === 'number' 
      ? width 
      : svgRef.current.clientWidth || 300;
    
    const containerHeight = height;
    
    // Dimensions of the chart area
    const chartWidth = containerWidth - margin.left - margin.right;
    const chartHeight = containerHeight - margin.top - margin.bottom;
    
    // Create scales for X and Y axes
    const xExtent = d3.extent(processedData, d => d.x) as [number, number];
    const yExtent = d3.extent(processedData, d => d.y) as [number, number];
    
    // Add a small padding to the extents to avoid points at the edges
    const xPadding = Math.max((xExtent[1] - xExtent[0]) * 0.05, 0.0001);
    const yPadding = Math.max((yExtent[1] - yExtent[0]) * 0.05, 0.0001);
    
    const xScale = d3.scaleLinear()
      .domain([xExtent[0] - xPadding, xExtent[1] + xPadding])
      .range([0, chartWidth])
      .nice();
    
    const yScale = d3.scaleLinear()
      .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
      .range([chartHeight, 0])
      .nice();
    
    // Create the main chart group
    const chartGroup = svg.append('g')
      .attr('class', 'chart-group')
      .attr('transform', `translate(${margin.left},${margin.top})`);
    
    // Add grid lines if enabled
    if (showGrid) {
      // Add X grid lines
      chartGroup.append('g')
        .attr('class', 'grid x-grid')
        .attr('transform', `translate(0,${chartHeight})`)
        .call(
          d3.axisBottom(xScale)
            .tickSize(-chartHeight)
            .tickFormat(() => '')
        )
        .selectAll('line')
        .attr('stroke', theme.palette.divider)
        .attr('stroke-opacity', 0.5);
      
      // Add Y grid lines
      chartGroup.append('g')
        .attr('class', 'grid y-grid')
        .call(
          d3.axisLeft(yScale)
            .tickSize(-chartWidth)
            .tickFormat(() => '')
        )
        .selectAll('line')
        .attr('stroke', theme.palette.divider)
        .attr('stroke-opacity', 0.5);
    }
    
    // Add X axis
    const xAxis = chartGroup.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${chartHeight})`)
      .call(d3.axisBottom(xScale).ticks(5).tickFormat(d => formatNumber(d as number)));
    
    // Add X axis label
    chartGroup.append('text')
      .attr('class', 'x-axis-label')
      .attr('text-anchor', 'middle')
      .attr('x', chartWidth / 2)
      .attr('y', chartHeight + margin.bottom - 5)
      .style('font-size', '0.75rem')
      .style('fill', theme.palette.text.secondary)
      .text(xAxisLabel || getPropertyDisplayName(xAxisProperty));
    
    // Add Y axis
    const yAxis = chartGroup.append('g')
      .attr('class', 'y-axis')
      .call(d3.axisLeft(yScale).ticks(5).tickFormat(d => formatNumber(d as number)));
    
    // Add Y axis label
    chartGroup.append('text')
      .attr('class', 'y-axis-label')
      .attr('text-anchor', 'middle')
      .attr('transform', 'rotate(-90)')
      .attr('x', -chartHeight / 2)
      .attr('y', -margin.left + 15)
      .style('font-size', '0.75rem')
      .style('fill', theme.palette.text.secondary)
      .text(yAxisLabel || getPropertyDisplayName(yAxisProperty));
    
    // Add scatter points
    const points = chartGroup.selectAll('.point')
      .data(processedData)
      .enter()
      .append('circle')
      .attr('class', 'point')
      .attr('cx', d => xScale(d.x))
      .attr('cy', d => yScale(d.y))
      .attr('r', 5)
      .attr('fill', color)
      .attr('opacity', 0.7)
      .attr('stroke', theme.palette.background.paper)
      .attr('stroke-width', 1)
      .style('cursor', onPointClick ? 'pointer' : 'default');
    
    // Add event handlers to points
    if (showTooltip) {
      points
        .on('mouseover', function(event, d) {
          showTooltipFunc(d, event);
        })
        .on('mouseout', hideTooltip);
    }
    
    if (onPointClick) {
      points.on('click', function(event, d) {
        onPointClick(d.originalData);
      });
    }
    
    // Set up zoom behavior if enabled
    if (allowZoom) {
      const zoom = d3.zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.5, 5])
        .extent([[0, 0], [containerWidth, containerHeight]])
        .on('zoom', handleZoom);
      
      svg.call(zoom)
        .call(
          // @ts-ignore (zoom.transform expects d3.ZoomBehavior which is hard to type correctly)
          zoom.transform,
          d3.zoomIdentity
            .translate(zoomState.translateX, zoomState.translateY)
            .scale(zoomState.scale)
        );
    }
    
    // Apply initial zoom state if not default
    if (allowZoom && (zoomState.scale !== 1 || zoomState.translateX !== 0 || zoomState.translateY !== 0)) {
      chartGroup.attr('transform', `translate(${zoomState.translateX + margin.left},${zoomState.translateY + margin.top}) scale(${zoomState.scale})`);
    }
    
    // Clean up function
    return () => {
      if (tooltipRef.current) {
        tooltipRef.current.style.display = 'none';
      }
    };
  }, [
    processedData, 
    windowSize, 
    color, 
    height, 
    width, 
    margin, 
    showGrid,
    xAxisLabel,
    yAxisLabel,
    zoomState,
    allowZoom,
    xAxisProperty,
    yAxisProperty,
    onPointClick,
    showTooltip
  ]);
  
  // Render skeleton if data is loading
  if (loading) {
    return (
      <Card className={className} style={style}>
        {title && <Typography variant="h6" gutterBottom>{title}</Typography>}
        <Skeleton variant="rectangular" width="100%" height={height} />
      </Card>
    );
  }
  
  return (
    <Card className={className} style={style}>
      {title && <Typography variant="h6" gutterBottom>{title}</Typography>}
      
      {/* Zoom controls */}
      {allowZoom && (
        <Box mb={1} display="flex" justifyContent="flex-end">
          <FormControl component="fieldset">
            <FormControlLabel
              control={
                <Switch
                  size="small"
                  checked={zoomState.scale !== 1 || zoomState.translateX !== 0 || zoomState.translateY !== 0}
                  onChange={resetZoom}
                />
              }
              label="Reset Zoom"
              labelPlacement="start"
            />
          </FormControl>
        </Box>
      )}
      
      {/* Chart container */}
      <Box width="100%" height={height} position="relative">
        <svg
          ref={svgRef}
          width={width}
          height={height}
          style={{ overflow: 'visible' }}
        />
        
        {/* Tooltip element (positioned absolutely) */}
        {showTooltip && (
          <div
            ref={tooltipRef}
            style={{
              position: 'fixed',
              display: 'none',
              backgroundColor: theme.palette.background.paper,
              border: `1px solid ${theme.palette.divider}`,
              borderRadius: 4,
              padding: 4,
              fontSize: '0.75rem',
              pointerEvents: 'none',
              zIndex: 1000,
              boxShadow: '0px 2px 8px rgba(0, 0, 0, 0.15)'
            }}
          />
        )}
      </Box>
    </Card>
  );
};

export default ScatterPlot;