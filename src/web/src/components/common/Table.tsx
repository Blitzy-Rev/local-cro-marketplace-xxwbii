import React, { useState, useEffect, useMemo, useCallback } from 'react'; // react 18.2+
import {
  Table as MuiTable,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableSortLabel,
  Paper,
  Checkbox,
  Box,
  Typography,
  CircularProgress
} from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { visuallyHidden } from '@mui/utils'; // v5.13+
import Pagination from './Pagination';
import theme from '../../theme';

// Type for sort direction
export type SortDirection = 'asc' | 'desc';

/**
 * Interface for table column configuration
 */
export interface TableColumn {
  /** Field name in the data object */
  field: string;
  /** Display name for the column header */
  headerName: string;
  /** Width of the column */
  width?: number | string;
  /** Text alignment within the column */
  align?: 'left' | 'center' | 'right';
  /** Whether this column can be sorted */
  sortable?: boolean;
  /** Custom render function for cell content */
  renderCell?: (row: any) => React.ReactNode;
  /** Custom render function for header content */
  renderHeader?: () => React.ReactNode;
  /** Function to extract value from row data */
  valueGetter?: (row: any) => any;
  /** Whether to hide this column */
  hide?: boolean;
  /** CSS class for cells in this column */
  cellClassName?: string | ((row: any) => string);
  /** CSS class for the header cell */
  headerClassName?: string;
}

/**
 * Props interface for table pagination
 */
export interface PaginationProps {
  /** Current page number (1-based) */
  page: number;
  /** Number of rows per page */
  pageSize: number;
  /** Total number of items across all pages */
  totalItems: number;
  /** Callback when page changes */
  onPageChange: (page: number) => void;
  /** Callback when page size changes */
  onPageSizeChange: (pageSize: number) => void;
  /** Available options for page size */
  pageSizeOptions?: number[];
  /** Whether to show page size selector */
  showPageSizeSelector?: boolean;
}

/**
 * Props interface for the Table component
 */
export interface TableProps {
  /** Array of data objects to display in the table */
  data: any[];
  /** Array of column definitions for the table */
  columns: TableColumn[];
  /** Pagination configuration for the table */
  pagination?: PaginationProps;
  /** Whether the table is in a loading state */
  loading?: boolean;
  /** Error message to display if data loading failed */
  error?: string | null;
  /** Whether rows can be selected with checkboxes */
  selectable?: boolean;
  /** Array of selected row IDs */
  selectedRows?: string[];
  /** Callback when row selection changes */
  onSelectionChange?: (selectedIds: string[]) => void;
  /** Field name to use as unique identifier for rows */
  idField?: string;
  /** Whether columns can be sorted */
  sortable?: boolean;
  /** Default field to sort by */
  defaultSortField?: string;
  /** Default sort direction */
  defaultSortDirection?: SortDirection;
  /** Callback when sort field or direction changes */
  onSortChange?: (field: string, direction: SortDirection) => void;
  /** Message to display when data is empty */
  emptyMessage?: string;
  /** Additional CSS class for the table container */
  className?: string;
  /** Additional inline styles for the table container */
  style?: React.CSSProperties;
  /** Whether the table header should stick to the top during scroll */
  stickyHeader?: boolean;
  /** Whether to use more compact row padding */
  dense?: boolean;
  /** Maximum height of the table with scrolling */
  maxHeight?: string | number;
}

// Styled components
const StyledTableContainer = styled(TableContainer)(({ theme }) => ({
  position: 'relative',
  overflow: 'auto',
  borderRadius: theme.shape.borderRadius,
  boxShadow: theme.shadows[1],
}));

const LoadingOverlay = styled(Box)({
  position: 'absolute',
  top: 0,
  left: 0,
  right: 0,
  bottom: 0,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  backgroundColor: 'rgba(255, 255, 255, 0.7)',
  zIndex: 1,
});

const ErrorContainer = styled(Box)(({ theme }) => ({
  padding: theme.spacing(2),
  textAlign: 'center',
  color: theme.palette.error.main,
}));

const EmptyMessageContainer = styled(Box)(({ theme }) => ({
  padding: theme.spacing(4),
  textAlign: 'center',
  color: theme.palette.text.secondary,
}));

/**
 * A reusable table component with sorting, pagination, and selection capabilities
 * 
 * @param props The component props
 * @returns The rendered table component
 */
const Table: React.FC<TableProps> = ({
  data = [],
  columns = [],
  pagination,
  loading = false,
  error = null,
  selectable = false,
  selectedRows = [],
  onSelectionChange,
  idField = 'id',
  sortable = false,
  defaultSortField = '',
  defaultSortDirection = 'asc',
  onSortChange,
  emptyMessage = 'No data available',
  className,
  style,
  stickyHeader = false,
  dense = false,
  maxHeight,
}) => {
  // State for sorting
  const [sortField, setSortField] = useState<string>(defaultSortField);
  const [sortDirection, setSortDirection] = useState<SortDirection>(defaultSortDirection);

  // State for selection
  const [selected, setSelected] = useState<string[]>(selectedRows || []);

  // Update selected state when selectedRows prop changes
  useEffect(() => {
    setSelected(selectedRows || []);
  }, [selectedRows]);

  // Handle sort change
  const handleSortChange = useCallback((field: string) => {
    const isAsc = sortField === field && sortDirection === 'asc';
    const newDirection = isAsc ? 'desc' : 'asc';
    
    setSortField(field);
    setSortDirection(newDirection);
    
    if (onSortChange) {
      onSortChange(field, newDirection);
    }
  }, [sortField, sortDirection, onSortChange]);

  // Check if a row is selected
  const isSelected = useCallback((id: string) => {
    return selected.indexOf(id) !== -1;
  }, [selected]);

  // Handle row selection
  const handleRowSelect = useCallback((id: string) => {
    const selectedIndex = selected.indexOf(id);
    let newSelected: string[] = [];

    if (selectedIndex === -1) {
      newSelected = [...selected, id];
    } else {
      newSelected = selected.filter(item => item !== id);
    }

    setSelected(newSelected);
    
    if (onSelectionChange) {
      onSelectionChange(newSelected);
    }
  }, [selected, onSelectionChange]);

  // Handle select all
  const handleSelectAll = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    if (event.target.checked) {
      const newSelected = data.map(row => row[idField]);
      setSelected(newSelected);
      if (onSelectionChange) {
        onSelectionChange(newSelected);
      }
      return;
    }
    
    setSelected([]);
    if (onSelectionChange) {
      onSelectionChange([]);
    }
  }, [data, idField, onSelectionChange]);

  // Visible columns (excluding hidden ones)
  const visibleColumns = useMemo(() => {
    return columns.filter(column => !column.hide);
  }, [columns]);

  // Compute number of selected items for optimized rendering
  const numSelected = useMemo(() => selected.length, [selected]);
  
  // Compute if all items are selected
  const allSelected = useMemo(() => {
    return data.length > 0 && numSelected === data.length;
  }, [data.length, numSelected]);
  
  // Compute if some items are selected
  const someSelected = useMemo(() => {
    return numSelected > 0 && numSelected < data.length;
  }, [numSelected, data.length]);

  return (
    <div className={className} style={style}>
      <StyledTableContainer 
        component={Paper} 
        sx={{ maxHeight: maxHeight }}
      >
        {loading && (
          <LoadingOverlay>
            <CircularProgress />
          </LoadingOverlay>
        )}
        
        <MuiTable
          stickyHeader={stickyHeader}
          size={dense ? 'small' : 'medium'}
          aria-label="data table"
        >
          <TableHead>
            <TableRow>
              {selectable && (
                <TableCell padding="checkbox">
                  <Checkbox
                    indeterminate={someSelected}
                    checked={allSelected}
                    onChange={handleSelectAll}
                    inputProps={{ 'aria-label': 'select all' }}
                    color="primary"
                  />
                </TableCell>
              )}
              
              {visibleColumns.map((column) => {
                const isSortable = sortable && column.sortable !== false;
                
                return (
                  <TableCell
                    key={column.field}
                    align={column.align || 'left'}
                    style={{ width: column.width }}
                    className={column.headerClassName}
                    sortDirection={sortField === column.field ? sortDirection : false}
                  >
                    {column.renderHeader ? (
                      column.renderHeader()
                    ) : isSortable ? (
                      <TableSortLabel
                        active={sortField === column.field}
                        direction={sortField === column.field ? sortDirection : 'asc'}
                        onClick={() => handleSortChange(column.field)}
                      >
                        {column.headerName}
                        {sortField === column.field ? (
                          <Box component="span" sx={visuallyHidden}>
                            {sortDirection === 'desc' ? 'sorted descending' : 'sorted ascending'}
                          </Box>
                        ) : null}
                      </TableSortLabel>
                    ) : (
                      column.headerName
                    )}
                  </TableCell>
                );
              })}
            </TableRow>
          </TableHead>
          
          {!error && (
            <TableBody>
              {data.length > 0 ? (
                data.map((row) => {
                  const rowId = row[idField];
                  const isItemSelected = isSelected(rowId);
                  
                  return (
                    <TableRow
                      hover
                      key={rowId}
                      selected={isItemSelected}
                      onClick={selectable ? () => handleRowSelect(rowId) : undefined}
                      role={selectable ? 'checkbox' : undefined}
                      aria-checked={selectable ? isItemSelected : undefined}
                      tabIndex={selectable ? -1 : undefined}
                    >
                      {selectable && (
                        <TableCell padding="checkbox">
                          <Checkbox
                            checked={isItemSelected}
                            color="primary"
                            onClick={(event) => event.stopPropagation()}
                            onChange={() => handleRowSelect(rowId)}
                          />
                        </TableCell>
                      )}
                      
                      {visibleColumns.map((column) => {
                        let cellValue;
                        
                        if (column.valueGetter) {
                          cellValue = column.valueGetter(row);
                        } else {
                          cellValue = row[column.field];
                        }

                        const cellClassName = typeof column.cellClassName === 'function'
                          ? column.cellClassName(row)
                          : column.cellClassName;
                        
                        return (
                          <TableCell
                            key={column.field}
                            align={column.align || 'left'}
                            className={cellClassName}
                          >
                            {column.renderCell ? column.renderCell(row) : cellValue}
                          </TableCell>
                        );
                      })}
                    </TableRow>
                  );
                })
              ) : (
                <TableRow>
                  <TableCell
                    colSpan={selectable ? visibleColumns.length + 1 : visibleColumns.length}
                    align="center"
                  >
                    <EmptyMessageContainer>
                      <Typography variant="body2">{emptyMessage}</Typography>
                    </EmptyMessageContainer>
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          )}
        </MuiTable>
        
        {error && (
          <ErrorContainer>
            <Typography variant="body2" color="error">
              {error}
            </Typography>
          </ErrorContainer>
        )}
      </StyledTableContainer>
      
      {pagination && (
        <Pagination
          page={pagination.page}
          pageSize={pagination.pageSize}
          totalItems={pagination.totalItems}
          onPageChange={pagination.onPageChange}
          onPageSizeChange={pagination.onPageSizeChange}
          pageSizeOptions={pagination.pageSizeOptions}
          showPageSizeSelector={pagination.showPageSizeSelector}
        />
      )}
    </div>
  );
};

export default Table;