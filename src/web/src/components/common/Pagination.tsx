import React, { useMemo } from 'react';
import { Box, Typography, IconButton, Select, MenuItem, FormControl, InputLabel } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { FirstPage, LastPage, NavigateBefore, NavigateNext } from '@mui/icons-material'; // v5.13+
import theme from '../../theme';

/**
 * Props interface for the Pagination component
 */
interface PaginationProps {
  /** Current page number (1-based) */
  page: number;
  /** Number of items per page */
  pageSize: number;
  /** Total number of items across all pages */
  totalItems: number;
  /** Callback when page changes */
  onPageChange: (page: number) => void;
  /** Callback when page size changes */
  onPageSizeChange: (pageSize: number) => void;
  /** Available options for page size selection */
  pageSizeOptions?: number[];
  /** Whether to show the page size selector */
  showPageSizeSelector?: boolean;
  /** Whether to show first/last page buttons */
  showFirstLastButtons?: boolean;
  /** Whether to show current range and total items */
  showPageInfo?: boolean;
  /** Whether the pagination controls are disabled */
  disabled?: boolean;
  /** Additional CSS class for the pagination container */
  className?: string;
}

// Styled components
const PaginationContainer = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'flex-end',
  padding: '16px 0',
  gap: '16px',
  flexWrap: 'wrap',
}));

const PageInfo = styled(Typography)(({ theme }) => ({
  color: theme.palette.text.secondary,
  fontSize: '0.875rem',
  whiteSpace: 'nowrap',
}));

const PageSizeSelector = styled(FormControl)(({ theme }) => ({
  minWidth: '120px',
  marginRight: '16px',
}));

const NavigationControls = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
}));

/**
 * A component that renders pagination controls for navigating through paginated data
 * 
 * @param props - The component props
 * @returns The rendered pagination component
 */
const Pagination: React.FC<PaginationProps> = ({
  page,
  pageSize,
  totalItems,
  onPageChange,
  onPageSizeChange,
  pageSizeOptions = [10, 25, 50, 100],
  showPageSizeSelector = true,
  showFirstLastButtons = true,
  showPageInfo = true,
  disabled = false,
  className,
}) => {
  // Calculate total pages
  const totalPages = useMemo(() => {
    return Math.max(1, Math.ceil(totalItems / pageSize));
  }, [totalItems, pageSize]);

  // Calculate current page range (e.g., "1-10 of 100")
  const pageRange = useMemo(() => {
    if (totalItems === 0) return '0-0 of 0';
    
    const startItem = (page - 1) * pageSize + 1;
    const endItem = Math.min(page * pageSize, totalItems);
    
    return `${startItem}-${endItem} of ${totalItems}`;
  }, [page, pageSize, totalItems]);

  // Handle navigation to first page
  const handleFirstPage = () => {
    if (page !== 1 && !disabled) {
      onPageChange(1);
    }
  };

  // Handle navigation to previous page
  const handlePreviousPage = () => {
    if (page > 1 && !disabled) {
      onPageChange(page - 1);
    }
  };

  // Handle navigation to next page
  const handleNextPage = () => {
    if (page < totalPages && !disabled) {
      onPageChange(page + 1);
    }
  };

  // Handle navigation to last page
  const handleLastPage = () => {
    if (page !== totalPages && !disabled) {
      onPageChange(totalPages);
    }
  };

  // Handle page size change
  const handlePageSizeChange = (event: React.ChangeEvent<{ value: unknown }>) => {
    if (!disabled) {
      const newPageSize = Number(event.target.value);
      onPageSizeChange(newPageSize);
      
      // Adjust current page to maintain approximate position in the data
      const currentFirstItem = (page - 1) * pageSize + 1;
      const newPage = Math.max(1, Math.ceil(currentFirstItem / newPageSize));
      if (newPage !== page) {
        onPageChange(newPage);
      }
    }
  };

  return (
    <PaginationContainer className={className}>
      {/* Page info section */}
      {showPageInfo && (
        <PageInfo variant="body2">
          {pageRange}
        </PageInfo>
      )}

      {/* Page size selector */}
      {showPageSizeSelector && (
        <PageSizeSelector size="small">
          <InputLabel id="pagination-page-size-label">Rows per page</InputLabel>
          <Select
            labelId="pagination-page-size-label"
            id="pagination-page-size"
            value={pageSize}
            onChange={handlePageSizeChange as any}
            label="Rows per page"
            disabled={disabled}
            size="small"
          >
            {pageSizeOptions.map((option) => (
              <MenuItem key={option} value={option}>
                {option}
              </MenuItem>
            ))}
          </Select>
        </PageSizeSelector>
      )}

      {/* Navigation buttons */}
      <NavigationControls>
        {/* First page button */}
        {showFirstLastButtons && (
          <IconButton
            onClick={handleFirstPage}
            disabled={page === 1 || disabled || totalItems === 0}
            aria-label="First page"
            size="small"
          >
            <FirstPage />
          </IconButton>
        )}

        {/* Previous page button */}
        <IconButton
          onClick={handlePreviousPage}
          disabled={page === 1 || disabled || totalItems === 0}
          aria-label="Previous page"
          size="small"
        >
          <NavigateBefore />
        </IconButton>

        {/* Page indicator */}
        <Typography variant="body2" sx={{ mx: 1 }}>
          {`${page} / ${totalPages}`}
        </Typography>

        {/* Next page button */}
        <IconButton
          onClick={handleNextPage}
          disabled={page >= totalPages || disabled || totalItems === 0}
          aria-label="Next page"
          size="small"
        >
          <NavigateNext />
        </IconButton>

        {/* Last page button */}
        {showFirstLastButtons && (
          <IconButton
            onClick={handleLastPage}
            disabled={page >= totalPages || disabled || totalItems === 0}
            aria-label="Last page"
            size="small"
          >
            <LastPage />
          </IconButton>
        )}
      </NavigationControls>
    </PaginationContainer>
  );
};

export default Pagination;