import React, { useState, useEffect, useMemo } from 'react';
import { Box, Grid, Typography, Divider, CircularProgress, Chip, Tooltip } from '@mui/material';
import { Science, Assessment } from '@mui/icons-material';

import { ResultData, ResultDetailed } from '../../../types/result';
import { Molecule } from '../../../types/molecule';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import Card from '../../../components/common/Card';
import useResults from '../hooks/useResults';
import { formatPropertyValue } from '../../../utils/formatters';

interface MoleculeResultsProps {
  result: ResultDetailed;
  className?: string;
}

interface MoleculeResultData {
  molecule?: Molecule;
  properties: { [key: string]: { value: string | number; unit?: string } };
}

const MoleculeResults: React.FC<MoleculeResultsProps> = ({ result, className }) => {
  const [loading, setLoading] = useState(true);
  const [moleculeData, setMoleculeData] = useState<Record<string, MoleculeResultData>>({});
  const { getResultData } = useResults();

  // Fetch result data for molecules when component mounts
  useEffect(() => {
    const fetchMoleculeData = async () => {
      if (!result || !result.id) {
        setLoading(false);
        return;
      }

      try {
        setLoading(true);
        const resultData = await getResultData(result.id);
        
        if (resultData) {
          // Group data points by molecule
          const groupedData: Record<string, MoleculeResultData> = {};
          
          resultData.forEach((dataPoint: ResultData) => {
            const { molecule_id, data_name, data_value, data_unit, molecule } = dataPoint;
            
            if (!groupedData[molecule_id]) {
              groupedData[molecule_id] = {
                molecule: molecule,
                properties: {}
              };
            }
            
            groupedData[molecule_id].properties[data_name] = {
              value: data_value,
              unit: data_unit
            };
          });
          
          setMoleculeData(groupedData);
        }
      } catch (error) {
        console.error('Error fetching molecule result data:', error);
      } finally {
        setLoading(false);
      }
    };

    fetchMoleculeData();
  }, [result, getResultData]);

  // Memo for molecules array to avoid re-rendering
  const molecules = useMemo(() => {
    return Object.keys(moleculeData).map(moleculeId => ({
      id: moleculeId,
      data: moleculeData[moleculeId]
    }));
  }, [moleculeData]);

  if (loading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: 200 }} className={className}>
        <CircularProgress />
        <Typography variant="body2" sx={{ ml: 2 }}>
          Loading molecule results...
        </Typography>
      </Box>
    );
  }

  if (molecules.length === 0) {
    return (
      <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'center', minHeight: 200, justifyContent: 'center' }} className={className}>
        <Assessment fontSize="large" color="action" />
        <Typography variant="body1" color="textSecondary" sx={{ mt: 2 }}>
          No molecule-specific result data available
        </Typography>
      </Box>
    );
  }

  return (
    <Box className={className}>
      <Box sx={{ mb: 3 }}>
        <Typography variant="h6" sx={{ display: 'flex', alignItems: 'center' }}>
          <Science sx={{ mr: 1 }} />
          Molecular Results Data
        </Typography>
        <Typography variant="body2" color="textSecondary">
          Showing experimental results for {molecules.length} molecule{molecules.length !== 1 ? 's' : ''}
        </Typography>
      </Box>
      
      <Grid container spacing={3}>
        {molecules.map(({ id, data }) => {
          const smiles = data.molecule?.smiles || '';
          const propertyEntries = Object.entries(data.properties);
          
          return (
            <Grid item xs={12} sm={6} md={4} key={id}>
              <Card>
                <Box sx={{ p: 2 }}>
                  <Typography variant="subtitle1" sx={{ mb: 1, fontWeight: 'bold' }}>
                    Molecule {id.slice(0, 8)}
                  </Typography>
                  
                  <Box sx={{ height: 200, mb: 2 }}>
                    <MoleculeViewer 
                      smiles={smiles}
                      molecule={data.molecule}
                      height={200}
                      width={300}
                      showControls={false}
                      showProperties={false}
                    />
                  </Box>
                  
                  <Divider sx={{ my: 2 }} />
                  
                  <Typography variant="subtitle2" sx={{ mb: 1 }}>
                    Experimental Results:
                  </Typography>
                  
                  {propertyEntries.length > 0 ? (
                    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
                      {propertyEntries.map(([name, { value, unit }]) => (
                        <Box 
                          key={name} 
                          sx={{ 
                            display: 'flex', 
                            justifyContent: 'space-between',
                            alignItems: 'center', 
                            py: 0.5
                          }}
                        >
                          <Tooltip title={name}>
                            <Typography variant="body2" sx={{ 
                              fontWeight: 'medium',
                              maxWidth: '60%',
                              overflow: 'hidden',
                              textOverflow: 'ellipsis',
                              whiteSpace: 'nowrap'
                            }}>
                              {name}:
                            </Typography>
                          </Tooltip>
                          <Chip 
                            label={formatPropertyValue(value, name, unit)}
                            size="small"
                            color="primary"
                            variant="outlined"
                          />
                        </Box>
                      ))}
                    </Box>
                  ) : (
                    <Typography variant="body2" color="textSecondary">
                      No property data available
                    </Typography>
                  )}
                </Box>
              </Card>
            </Grid>
          );
        })}
      </Grid>
    </Box>
  );
};

export default MoleculeResults;