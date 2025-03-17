import React, { useState, useEffect, useRef, useMemo } from 'react'; // react version 18.2+
import * as d3 from 'd3'; // d3 version 7.8+
import { 
  Box, 
  Typography, 
  Tooltip, 
  Skeleton,
  FormControl,
  FormLabel,
  Select,
  MenuItem,
  Switch,
  FormControlLabel
} from '@mui/material'; // @mui/material version 5.13+

// Internal imports
import useWindowSize from '../../hooks/useWindowSize';
import { formatNumber, formatPropertyValue } from '../../utils/formatters';
import { getPropertyDisplayName, getPropertyUnit } from '../../utils/molecularUtils';
import theme from '../../theme';
import Card from '../common/Card';

/**
 * Props interface for the PropertyChart component
 */
interface PropertyChartProps {
  /** Array of data objects containing properties to chart */
  data: Array<{ [key: string]: number | string }>;
  /** Chart title */
  title?: string;
  /** X-axis label */
  xAxisLabel?: string;
  /** Y-axis label */
  yAxisLabel?: string;
  /** Property name to visualize */
  propertyName: string;
  /** Property to use for x-axis values (default is index) */
  xProperty?: string;
  /** Chart type: 'line', 'bar', or 'area' */
  chartType?: string;
  /** Array of colors for multi-series charts */
  colors?: string[];
  /** Chart height in pixels */
  height?: number;
  /** Chart width (number or string like '100%') */
  width?: number | string;
  /** Chart margins */
  margin?: { top: number; right: number; bottom: number; left: number };
  /** Whether to show tooltips on hover */
  showTooltip?: boolean;
  /** Whether to show legend */
  showLegend?: boolean;
  /** Whether to show grid lines */
  showGrid?: boolean;
  /** Whether chart is in loading state */
  loading?: boolean;
  /** Whether to allow changing chart type */
  allowChartTypeChange?: boolean;
  /** Handler for point click event */
  onPointClick?: (point: any) => void;
  /** Additional class name */
  className?: string;
  /** Additional styles */
  style?: React.CSSProperties;
}

/**
 * Interface representing a data point in the chart
 */
interface ChartPoint {
  /** X-coordinate value */
  x: number;
  /** Y-coordinate value */
  y: number;
  /** Label for the data point */
  label: string;
  /** Original data object */
  originalData: any;
}

/**
 * A component that renders various chart types for visualizing molecular property data
 */
const PropertyChart: React.FC<PropertyChartProps> = ({
  data,
  title,
  xAxisLabel,
  yAxisLabel,
  propertyName,
  xProperty,
  chartType = 'line',
  colors = [
    theme.palette.primary.main,
    theme.palette.secondary.main, 
    theme.palette.success.main,
    theme.palette.error.main
  ],
  height = 300,
  width = '100%',
  margin = { top: 20, right: 20, bottom: 50, left: 50 },
  showTooltip = true,
  showLegend = true,
  showGrid = true,
  loading = false,
  allowChartTypeChange = false,
  onPointClick,
  className,
  style
}) => {
  // Refs for D3 manipulation
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  
  // Get window size for responsive resizing
  const windowSize = useWindowSize();
  
  // State for current chart type and options
  const [currentChartType, setCurrentChartType] = useState<string>(chartType);
  const [showGridLines, setShowGridLines] = useState<boolean>(showGrid);
  const [showTooltips, setShowTooltips] = useState<boolean>(showTooltip);
  
  // Process the data into format needed for visualization
  const processedData = useMemo(() => {
    if (!data || data.length === 0) return [];
    
    return data.map((item, index) => {
      const x = xProperty ? Number(item[xProperty]) : index;
      const y = Number(item[propertyName]);
      const label = xProperty ? String(item[xProperty]) : String(index + 1);
      
      return {
        x,
        y: isNaN(y) ? 0 : y,
        label,
        originalData: item
      };
    }).filter(point => !isNaN(point.x) && !isNaN(point.y));
  }, [data, propertyName, xProperty]);
  
  // Handle chart type change
  const handleChartTypeChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setCurrentChartType(event.target.value);
  };
  
  // Handle grid lines toggle
  const handleGridToggle = (event: React.ChangeEvent<HTMLInputElement>) => {
    setShowGridLines(event.target.checked);
  };
  
  // Handle tooltip toggle
  const handleTooltipToggle = (event: React.ChangeEvent<HTMLInputElement>) => {
    setShowTooltips(event.target.checked);
  };
  
  // Handle click on data points
  const handlePointClick = (point: ChartPoint, event: React.MouseEvent) => {
    if (onPointClick) {
      onPointClick(point.originalData);
    }
  };
  
  // Show tooltip on hover
  const showTooltipHandler = (point: ChartPoint, event: React.MouseEvent) => {
    if (!showTooltips || !tooltipRef.current) return;
    
    const tooltip = tooltipRef.current;
    tooltip.style.opacity = '1';
    tooltip.style.left = `${event.clientX + 10}px`;
    tooltip.style.top = `${event.clientY - 10}px`;
    
    const xValue = xProperty ? point.label : `Item ${point.label}`;
    const yValue = formatPropertyValue(propertyName, point.y);
    const unit = getPropertyUnit(propertyName);
    
    tooltip.innerHTML = `
      <div><strong>${xValue}</strong></div>
      <div>${getPropertyDisplayName(propertyName)}: ${yValue}${unit ? ` ${unit}` : ''}</div>
    `;
  };
  
  // Hide tooltip
  const hideTooltipHandler = () => {
    if (tooltipRef.current) {
      tooltipRef.current.style.opacity = '0';
    }
  };
  
  // Create and update the chart
  useEffect(() => {
    // Skip if no data or no svg element
    if (!svgRef.current || processedData.length === 0) return;
    
    // Clear existing content
    d3.select(svgRef.current).selectAll('*').remove();
    
    // Get container dimensions
    const containerWidth = typeof width === 'number' 
      ? width 
      : svgRef.current.parentElement?.clientWidth || 300;
    
    // Calculate chart dimensions
    const chartWidth = containerWidth - margin.left - margin.right;
    const chartHeight = height - margin.top - margin.bottom;
    
    // Create SVG
    const svg = d3.select(svgRef.current)
      .attr('width', containerWidth)
      .attr('height', height);
    
    // Create chart group with margin
    const chart = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);
    
    // Create scales
    const xMin = d3.min(processedData, d => d.x) as number;
    const xMax = d3.max(processedData, d => d.x) as number;
    
    const xScale = d3.scaleLinear()
      .domain([
        xMin > 0 && xProperty ? 0 : xMin, 
        xMax
      ])
      .range([0, chartWidth])
      .nice();
    
    const yMin = 0;
    const yMax = d3.max(processedData, d => d.y) as number;
    
    const yScale = d3.scaleLinear()
      .domain([yMin, yMax * 1.05]) // Add 5% padding at the top
      .range([chartHeight, 0])
      .nice();
    
    // Create clip path to prevent drawing outside chart area
    chart.append('clipPath')
      .attr('id', 'chart-area-clip')
      .append('rect')
      .attr('width', chartWidth)
      .attr('height', chartHeight)
      .attr('x', 0)
      .attr('y', 0);
    
    // Create axes
    const xAxis = d3.axisBottom(xScale)
      .ticks(5)
      .tickFormat(d => formatNumber(d as number));
    
    const yAxis = d3.axisLeft(yScale)
      .ticks(5)
      .tickFormat(d => formatNumber(d as number));
    
    // Add x-axis
    chart.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${chartHeight})`)
      .call(xAxis)
      .selectAll('text')
      .attr('font-size', '12px');
    
    // Add y-axis
    chart.append('g')
      .attr('class', 'y-axis')
      .call(yAxis)
      .selectAll('text')
      .attr('font-size', '12px');
    
    // Add axis labels
    if (xAxisLabel) {
      chart.append('text')
        .attr('class', 'x-axis-label')
        .attr('text-anchor', 'middle')
        .attr('x', chartWidth / 2)
        .attr('y', chartHeight + 40)
        .attr('font-size', '14px')
        .text(xAxisLabel);
    }
    
    if (yAxisLabel) {
      chart.append('text')
        .attr('class', 'y-axis-label')
        .attr('text-anchor', 'middle')
        .attr('transform', 'rotate(-90)')
        .attr('x', -chartHeight / 2)
        .attr('y', -40)
        .attr('font-size', '14px')
        .text(yAxisLabel || getPropertyDisplayName(propertyName));
    }
    
    // Add grid lines if enabled
    if (showGridLines) {
      // Add x grid lines
      chart.append('g')
        .attr('class', 'grid x-grid')
        .attr('transform', `translate(0,${chartHeight})`)
        .call(
          d3.axisBottom(xScale)
            .tickSize(-chartHeight)
            .tickFormat(() => '')
        )
        .selectAll('line')
        .attr('stroke', 'rgba(0,0,0,0.1)');
      
      // Add y grid lines
      chart.append('g')
        .attr('class', 'grid y-grid')
        .call(
          d3.axisLeft(yScale)
            .tickSize(-chartWidth)
            .tickFormat(() => '')
        )
        .selectAll('line')
        .attr('stroke', 'rgba(0,0,0,0.1)');
      
      // Hide tick lines as they're redundant with the grid
      chart.selectAll('.grid .domain').style('display', 'none');
    }
    
    // Create chart content based on type
    const chartContent = chart.append('g')
      .attr('clip-path', 'url(#chart-area-clip)');
    
    // Render the appropriate chart type
    switch (currentChartType) {
      case 'bar':
        renderBarChart(chartContent, xScale, yScale);
        break;
      case 'area':
        renderAreaChart(chartContent, xScale, yScale);
        break;
      case 'line':
      default:
        renderLineChart(chartContent, xScale, yScale);
        break;
    }
    
    // Add legend if enabled
    if (showLegend) {
      renderLegend(chart);
    }
    
    // Helper function to render line chart
    function renderLineChart(
      chartGroup: d3.Selection<SVGGElement, unknown, null, undefined>, 
      xScale: d3.ScaleLinear<number, number>, 
      yScale: d3.ScaleLinear<number, number>
    ) {
      // Sort data by x value for proper line rendering
      const sortedData = [...processedData].sort((a, b) => a.x - b.x);
      
      // Create line generator
      const line = d3.line<ChartPoint>()
        .x(d => xScale(d.x))
        .y(d => yScale(d.y))
        .curve(d3.curveMonotoneX);
      
      // Add line path
      chartGroup.append('path')
        .datum(sortedData)
        .attr('class', 'line')
        .attr('fill', 'none')
        .attr('stroke', colors[0])
        .attr('stroke-width', 2)
        .attr('d', line);
      
      // Add data points
      chartGroup.selectAll('.data-point')
        .data(processedData)
        .enter()
        .append('circle')
        .attr('class', 'data-point')
        .attr('cx', d => xScale(d.x))
        .attr('cy', d => yScale(d.y))
        .attr('r', 5)
        .attr('fill', colors[0])
        .style('cursor', onPointClick ? 'pointer' : 'default')
        .on('mouseover', function(event, d) {
          d3.select(this).attr('r', 7);
          showTooltipHandler(d, event);
        })
        .on('mouseout', function() {
          d3.select(this).attr('r', 5);
          hideTooltipHandler();
        })
        .on('click', function(event, d) {
          handlePointClick(d, event);
        });
    }
    
    // Helper function to render bar chart
    function renderBarChart(
      chartGroup: d3.Selection<SVGGElement, unknown, null, undefined>, 
      xScale: d3.ScaleLinear<number, number>, 
      yScale: d3.ScaleLinear<number, number>
    ) {
      // Calculate bar width based on data size and chart width
      const barPadding = 0.2; // 20% padding between bars
      const barWidth = Math.min(
        chartWidth / processedData.length * (1 - barPadding),
        50 // Maximum bar width
      );
      
      // Add bars
      chartGroup.selectAll('.bar')
        .data(processedData)
        .enter()
        .append('rect')
        .attr('class', 'bar')
        .attr('x', d => xScale(d.x) - barWidth / 2)
        .attr('y', d => yScale(d.y))
        .attr('width', barWidth)
        .attr('height', d => chartHeight - yScale(d.y))
        .attr('fill', colors[0])
        .attr('rx', 2) // Rounded corners
        .style('cursor', onPointClick ? 'pointer' : 'default')
        .on('mouseover', function(event, d) {
          d3.select(this)
            .attr('fill', d3.color(colors[0])?.darker(0.2) as string);
          showTooltipHandler(d, event);
        })
        .on('mouseout', function() {
          d3.select(this).attr('fill', colors[0]);
          hideTooltipHandler();
        })
        .on('click', function(event, d) {
          handlePointClick(d, event);
        });
    }
    
    // Helper function to render area chart
    function renderAreaChart(
      chartGroup: d3.Selection<SVGGElement, unknown, null, undefined>, 
      xScale: d3.ScaleLinear<number, number>, 
      yScale: d3.ScaleLinear<number, number>
    ) {
      // Sort data by x value for proper area rendering
      const sortedData = [...processedData].sort((a, b) => a.x - b.x);
      
      // Create area generator
      const area = d3.area<ChartPoint>()
        .x(d => xScale(d.x))
        .y0(chartHeight)
        .y1(d => yScale(d.y))
        .curve(d3.curveMonotoneX);
      
      // Add area path
      chartGroup.append('path')
        .datum(sortedData)
        .attr('class', 'area')
        .attr('fill', d3.color(colors[0])?.copy({opacity: 0.5}) as string)
        .attr('d', area);
      
      // Add line on top of area
      const line = d3.line<ChartPoint>()
        .x(d => xScale(d.x))
        .y(d => yScale(d.y))
        .curve(d3.curveMonotoneX);
      
      chartGroup.append('path')
        .datum(sortedData)
        .attr('class', 'line')
        .attr('fill', 'none')
        .attr('stroke', colors[0])
        .attr('stroke-width', 2)
        .attr('d', line);
      
      // Add data points
      chartGroup.selectAll('.data-point')
        .data(processedData)
        .enter()
        .append('circle')
        .attr('class', 'data-point')
        .attr('cx', d => xScale(d.x))
        .attr('cy', d => yScale(d.y))
        .attr('r', 5)
        .attr('fill', colors[0])
        .style('cursor', onPointClick ? 'pointer' : 'default')
        .on('mouseover', function(event, d) {
          d3.select(this).attr('r', 7);
          showTooltipHandler(d, event);
        })
        .on('mouseout', function() {
          d3.select(this).attr('r', 5);
          hideTooltipHandler();
        })
        .on('click', function(event, d) {
          handlePointClick(d, event);
        });
    }
    
    // Helper function to render legend
    function renderLegend(
      chartGroup: d3.Selection<SVGGElement, unknown, null, undefined>
    ) {
      const legend = chartGroup.append('g')
        .attr('class', 'legend')
        .attr('transform', `translate(${chartWidth - 120}, 10)`);
      
      // Add background rectangle for the legend
      legend.append('rect')
        .attr('width', 120)
        .attr('height', 30)
        .attr('fill', 'white')
        .attr('opacity', 0.8)
        .attr('rx', 4); // Rounded corners
      
      // Add color swatch
      legend.append('rect')
        .attr('x', 10)
        .attr('y', 10)
        .attr('width', 12)
        .attr('height', 12)
        .attr('fill', colors[0]);
      
      // Add legend text
      legend.append('text')
        .attr('x', 30)
        .attr('y', 20)
        .attr('font-size', '12px')
        .text(getPropertyDisplayName(propertyName));
    }
  }, [
    processedData, 
    currentChartType, 
    showGridLines, 
    showTooltips, 
    windowSize, 
    colors, 
    height, 
    width, 
    margin, 
    xAxisLabel, 
    yAxisLabel, 
    showLegend, 
    propertyName, 
    onPointClick, 
    xProperty
  ]);
  
  // Render loading skeleton if data is loading
  if (loading) {
    return (
      <Card className={className} style={style}>
        {title && (
          <Box p={2} pb={0}>
            <Typography variant="h6">{title}</Typography>
          </Box>
        )}
        <Box p={2}>
          <Skeleton variant="rectangular" width="100%" height={height} />
        </Box>
      </Card>
    );
  }
  
  // Render empty state if no data
  if (!processedData.length) {
    return (
      <Card className={className} style={style}>
        {title && (
          <Box p={2} pb={0}>
            <Typography variant="h6">{title}</Typography>
          </Box>
        )}
        <Box 
          p={2} 
          display="flex" 
          alignItems="center" 
          justifyContent="center" 
          height={height}
          bgcolor="background.default"
          borderRadius={1}
        >
          <Typography color="text.secondary">
            No data available for this property
          </Typography>
        </Box>
      </Card>
    );
  }
  
  return (
    <Card className={className} style={style}>
      <Box p={2} pb={allowChartTypeChange ? 1 : 0}>
        <Box display="flex" justifyContent="space-between" alignItems="center">
          {title && <Typography variant="h6">{title}</Typography>}
          
          {allowChartTypeChange && (
            <FormControl size="small" variant="outlined">
              <Select
                value={currentChartType}
                onChange={handleChartTypeChange}
                displayEmpty
                inputProps={{ 'aria-label': 'chart type' }}
              >
                <MenuItem value="line">Line Chart</MenuItem>
                <MenuItem value="bar">Bar Chart</MenuItem>
                <MenuItem value="area">Area Chart</MenuItem>
              </Select>
            </FormControl>
          )}
        </Box>
        
        {allowChartTypeChange && (
          <Box display="flex" justifyContent="flex-end" mt={1}>
            <FormControlLabel
              control={
                <Switch
                  checked={showGridLines}
                  size="small"
                  onChange={handleGridToggle}
                />
              }
              label="Grid"
            />
            <FormControlLabel
              control={
                <Switch
                  checked={showTooltips}
                  size="small"
                  onChange={handleTooltipToggle}
                />
              }
              label="Tooltip"
            />
          </Box>
        )}
      </Box>
      
      <Box p={2} pt={0} position="relative">
        <svg 
          ref={svgRef} 
          width="100%" 
          height={height}
          style={{ overflow: 'visible' }}
        ></svg>
        <div
          ref={tooltipRef}
          style={{
            position: 'fixed',
            opacity: 0,
            pointerEvents: 'none',
            backgroundColor: 'white',
            padding: '8px',
            border: '1px solid #ccc',
            borderRadius: '4px',
            boxShadow: '0 2px 5px rgba(0,0,0,0.15)',
            zIndex: 1000,
            transition: 'opacity 0.2s',
            fontSize: '14px',
            maxWidth: '200px'
          }}
        ></div>
      </Box>
    </Card>
  );
};

export default PropertyChart;