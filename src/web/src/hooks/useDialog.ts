import { useState, useCallback } from 'react'; // React 18.2+

/**
 * A custom React hook that provides state management and handlers for controlling
 * dialog visibility throughout the application. It simplifies the implementation
 * of modal dialogs by abstracting the open/close logic.
 * 
 * @param initialState - Optional boolean indicating whether the dialog should start open (default: false)
 * @returns An object containing the dialog state and control functions
 * 
 * @example
 * // Basic usage
 * const { open, handleOpen, handleClose } = useDialog();
 * 
 * // With an initially open dialog
 * const { open, handleOpen, handleClose, toggleOpen } = useDialog(true);
 * 
 * // In a component
 * return (
 *   <>
 *     <Button onClick={handleOpen}>Open Dialog</Button>
 *     <Dialog open={open} onClose={handleClose}>
 *       {/* Dialog content */}
 *     </Dialog>
 *   </>
 * );
 */
export const useDialog = (initialState = false) => {
  // State to track dialog open/closed status
  const [open, setOpen] = useState<boolean>(initialState);

  // Handler to open the dialog
  const handleOpen = useCallback(() => {
    setOpen(true);
  }, []);

  // Handler to close the dialog
  const handleClose = useCallback(() => {
    setOpen(false);
  }, []);

  // Handler to toggle the dialog's current state
  const toggleOpen = useCallback(() => {
    setOpen(prevOpen => !prevOpen);
  }, []);

  // Return an object with the state and all control functions
  return {
    open,
    handleOpen,
    handleClose,
    toggleOpen
  };
};

export default useDialog;