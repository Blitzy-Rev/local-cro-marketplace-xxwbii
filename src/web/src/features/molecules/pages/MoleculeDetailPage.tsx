import React, { useState, useEffect } from 'react'; // React 18.2+
import { useParams, useNavigate } from 'react-router-dom'; // react-router-dom v6.10+
import { Box, Container, Typography, Paper } from '@mui/material'; // @mui/material v5.13+
import { styled } from '@mui/material/styles'; // @mui/material v5.13+

import MoleculeDetail from '../components/MoleculeDetail';
import useMolecules from '../hooks/useMolecules';
import Loading from '../../../components/common/Loading';
import useToast from '../../../hooks/useToast';
import { getMoleculeById } from '../../../api/molecules';
import { Molecule, MoleculeDetailResponse, Library, Experiment, FlagStatus } from '../../../types';

/**
 * Interface for the API response containing molecule details
 */
interface MoleculeDetailResponse {
  molecule: Molecule;
  libraries: Library[];
  experiments: Experiment[];
}

/**
 * Styled container for the MoleculeDetailPage to manage layout
 */
const StyledContainer = styled(Container)({
  maxWidth: 'lg',
  marginTop: theme.spacing(4),
  marginBottom: theme.spacing(4),
});

/**
 * Styled paper component for consistent styling
 */
const DetailPaper = styled(Paper)({
  padding: theme.spacing(3),
  marginBottom: theme.spacing(3),
});

/**
 * Page component that displays detailed information about a specific molecule
 */
const MoleculeDetailPage: React.FC = () => {
  // LD1: Extract molecule ID from URL parameters using useParams hook
  const { id } = useParams<{ id: string }>();

  // LD1: Initialize state for molecule data, libraries, experiments, and loading status
  const [molecule, setMolecule] = useState<Molecule | null>(null);
  const [libraries, setLibraries] = useState<Library[]>([]);
  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);

  // LD1: Initialize navigation function using useNavigate hook
  const navigate = useNavigate();

  // LD1: Initialize toast notification hook using useToast
  const { showToast } = useToast();

  // LD1: Initialize molecule management functions using useMolecules hook
  const { updateMoleculeFlag, addToLibrary, addToExperiment } = useMolecules();

  // LD1: Fetch molecule details when component mounts or ID changes
  useEffect(() => {
    const fetchMoleculeDetails = async () => {
      setLoading(true);
      setError(null);
      try {
        if (id) {
          const response = await getMoleculeById(id);
          if (response.data) {
            const moleculeData = response.data.data as MoleculeDetailResponse;
            setMolecule(moleculeData.molecule);
            setLibraries(moleculeData.libraries);
            setExperiments(moleculeData.experiments);
          } else {
            setError('Molecule not found');
          }
        }
      } catch (e: any) {
        setError(e.message || 'Failed to fetch molecule details');
      } finally {
        setLoading(false);
      }
    };

    if (id) {
      fetchMoleculeDetails();
    }
  }, [id]);

  // LD1: Handle loading state by displaying a loading indicator
  if (loading) {
    return (
      <StyledContainer>
        <Loading message="Loading molecule details..." fullScreen />
      </StyledContainer>
    );
  }

  // LD1: Handle error state by displaying an error message
  if (error) {
    return (
      <StyledContainer>
        <DetailPaper>
          <Typography color="error">{error}</Typography>
        </DetailPaper>
      </StyledContainer>
    );
  }

  // LD1: Implement handleBack function to navigate back to molecules list
  const handleBack = () => {
    navigate('/molecules');
  };

  // LD1: Implement handleFlagChange function to update molecule flag status
  const handleFlagChange = async (moleculeToUpdate: Molecule, newFlagStatus: FlagStatus | null) => {
    if (moleculeToUpdate && moleculeToUpdate.id) {
      await updateMoleculeFlag(moleculeToUpdate.id, newFlagStatus || '');
    } else {
      showToast({
        type: 'error',
        message: 'Failed to update molecule flag: Invalid molecule data',
      });
    }
  };

  // LD1: Implement handleAddToLibrary function to add molecule to a library
  const handleAddToLibrary = async (moleculeToAdd: Molecule) => {
    if (moleculeToAdd) {
      // You might want to prompt the user to select a library here
      // For now, we'll just show a toast message
      showToast({
        type: 'info',
        message: `Adding molecule ${moleculeToAdd.smiles} to a library...`,
      });
      // In a real implementation, you would dispatch an action to open a modal
      // where the user can select the library to add the molecule to.
    }
  };

  // LD1: Implement handleAddToExperiment function to add molecule to an experiment
  const handleAddToExperiment = async (moleculeToAdd: Molecule) => {
    if (moleculeToAdd) {
      // You might want to prompt the user to select an experiment here
      // For now, we'll just show a toast message
      showToast({
        type: 'info',
        message: `Adding molecule ${moleculeToAdd.smiles} to an experiment...`,
      });
      // In a real implementation, you would dispatch an action to open a modal
      // where the user can select the experiment to add the molecule to.
    }
  };

  // LD1: Render MoleculeDetail component with molecule data and handler functions
  return (
    <StyledContainer>
      {molecule && (
        <MoleculeDetail
          molecule={molecule}
          libraries={libraries}
          experiments={experiments}
          onBack={handleBack}
          onFlagChange={handleFlagChange}
          onAddToLibrary={handleAddToLibrary}
          onAddToExperiment={handleAddToExperiment}
        />
      )}
    </StyledContainer>
  );
};

export default MoleculeDetailPage;