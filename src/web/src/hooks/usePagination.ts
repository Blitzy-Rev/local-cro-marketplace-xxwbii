import { useState, useEffect, useCallback, useMemo } from 'react'; // v18.2+

/**
 * Options for the usePagination hook
 */
interface UsePaginationOptions {
  /**
   * The initial page number (defaults to 1)
   */
  initialPage?: number;
  
  /**
   * The initial number of items per page (defaults to 10)
   */
  initialPageSize?: number;
  
  /**
   * The total number of items in the dataset (defaults to 0)
   */
  totalItems?: number;
}

/**
 * Object containing information about the current page range
 */
interface PageInfo {
  /**
   * The index of the first item on the current page
   */
  startItem: number;
  
  /**
   * The index of the last item on the current page
   */
  endItem: number;
  
  /**
   * The total number of items in the dataset
   */
  totalItems: number;
}

/**
 * Result object containing pagination state and control functions
 */
interface PaginationResult {
  /**
   * The current page number
   */
  page: number;
  
  /**
   * The current number of items per page
   */
  pageSize: number;
  
  /**
   * The total number of pages based on total items and page size
   */
  totalPages: number;
  
  /**
   * Boolean indicating if there is a next page
   */
  hasNextPage: boolean;
  
  /**
   * Boolean indicating if there is a previous page
   */
  hasPreviousPage: boolean;
  
  /**
   * Function to navigate to a specific page
   * @param page The page number to navigate to
   */
  goToPage: (page: number) => void;
  
  /**
   * Function to navigate to the next page
   */
  nextPage: () => void;
  
  /**
   * Function to navigate to the previous page
   */
  previousPage: () => void;
  
  /**
   * Function to navigate to the first page
   */
  firstPage: () => void;
  
  /**
   * Function to navigate to the last page
   */
  lastPage: () => void;
  
  /**
   * Function to change the number of items per page
   * @param size The new page size
   */
  setPageSize: (size: number) => void;
  
  /**
   * Object containing information about the current page range
   */
  pageInfo: PageInfo;
}

/**
 * A hook that provides pagination state and controls for navigating through paginated data.
 * This hook is designed to work with API endpoints that support pagination and is used
 * throughout the application for handling large datasets.
 * 
 * @param options - Configuration options for the pagination
 * @returns An object containing pagination state and control functions
 * 
 * @example
 * ```tsx
 * const {
 *   page,
 *   pageSize,
 *   totalPages,
 *   hasNextPage,
 *   hasPreviousPage,
 *   goToPage,
 *   nextPage,
 *   previousPage,
 *   firstPage,
 *   lastPage,
 *   setPageSize,
 *   pageInfo
 * } = usePagination({
 *   initialPage: 1,
 *   initialPageSize: 25,
 *   totalItems: 100
 * });
 * ```
 */
export default function usePagination({
  initialPage = 1,
  initialPageSize = 10,
  totalItems = 0,
}: UsePaginationOptions = {}): PaginationResult {
  // Initialize state for current page and page size with provided initial values
  const [page, setPage] = useState(initialPage);
  const [pageSize, setPageSizeState] = useState(initialPageSize);

  // Calculate total pages based on total items and page size
  const totalPages = useMemo(() => {
    return totalItems > 0 ? Math.max(1, Math.ceil(totalItems / pageSize)) : 1;
  }, [totalItems, pageSize]);

  // Calculate if there are next and previous pages based on current page and total pages
  const hasNextPage = page < totalPages;
  const hasPreviousPage = page > 1;

  // Create memoized page info object with current range and total items
  const pageInfo = useMemo<PageInfo>(() => {
    const startItem = totalItems > 0 ? (page - 1) * pageSize + 1 : 0;
    const endItem = Math.min(startItem + pageSize - 1, totalItems);
    return {
      startItem,
      endItem,
      totalItems,
    };
  }, [page, pageSize, totalItems]);

  // Implement goToPage function to navigate to a specific page with validation
  const goToPage = useCallback((newPage: number) => {
    const validatedPage = Math.max(1, Math.min(newPage, totalPages));
    setPage(validatedPage);
  }, [totalPages]);

  // Implement nextPage function to increment page if hasNextPage is true
  const nextPage = useCallback(() => {
    if (hasNextPage) {
      setPage((prevPage) => prevPage + 1);
    }
  }, [hasNextPage]);

  // Implement previousPage function to decrement page if hasPreviousPage is true
  const previousPage = useCallback(() => {
    if (hasPreviousPage) {
      setPage((prevPage) => prevPage - 1);
    }
  }, [hasPreviousPage]);

  // Implement firstPage function to navigate to page 1
  const firstPage = useCallback(() => {
    setPage(1);
  }, []);

  // Implement lastPage function to navigate to the last page
  const lastPage = useCallback(() => {
    setPage(totalPages);
  }, [totalPages]);

  // Implement setPageSize function to update page size and recalculate current page
  const setPageSize = useCallback((newSize: number) => {
    if (newSize < 1) return;
    
    // Adjust current page to maintain the first visible item when possible
    const firstItemIndex = (page - 1) * pageSize;
    const newPage = Math.floor(firstItemIndex / newSize) + 1;
    
    setPageSizeState(newSize);
    setPage(Math.min(newPage, Math.ceil(totalItems / newSize) || 1));
  }, [page, pageSize, totalItems]);

  // Ensure current page is valid when total items or total pages change
  useEffect(() => {
    if (page > totalPages && totalPages > 0) {
      setPage(totalPages);
    }
  }, [totalPages, page]);

  // Return an object containing all pagination state and control functions
  return {
    page,
    pageSize,
    totalPages,
    hasNextPage,
    hasPreviousPage,
    goToPage,
    nextPage,
    previousPage,
    firstPage,
    lastPage,
    setPageSize,
    pageInfo,
  };
}