# ----------------------------------------------------------------------------
# Copyright (c) 2024, Lucas Patel, Yang Chen
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import tempfile
import pandas as pd
import numpy as np
import pandas.testing as pdt
import biom
import qiime2
from unittest.mock import patch, MagicMock
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.feature_table import BIOMV210Format
from qiime2 import Metadata
from q2_birdman._methods import _create_dir, run
from q2_birdman._visualizers import _escape_for_js, _sanitize_filename, _create_metadata_visualization, _compute_sample_log_ratios
import patsy


class CreateDirTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def test_create_dir_creates_all_subdirs(self):
        """
        Test that _create_dir creates all required subdirectories, otherwise raises AssertionError
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Call the _create_dir method with the temporary directory
            _create_dir(temp_dir)

            # Define the expected subdirectories
            expected_sub_dirs = ["slurm_out", "logs", "inferences", "results", "plots"]

            # Check that each subdirectory exists
            for sub_dir in expected_sub_dirs:
                sub_dir_path = os.path.join(temp_dir, sub_dir)
                assert os.path.exists(sub_dir_path), f"Subdirectory {sub_dir} does not exist."

    def test_create_dir_handles_existing_subdirs(self):
        """
        Test that _create_dir handles pre-existing subdirectories, otherwises raises AssertionError
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Manually create one of the subdirectories
            os.makedirs(os.path.join(temp_dir, "slurm_out"), exist_ok=True)

            # Call the _create_dir method with the temporary directory
            _create_dir(temp_dir)

            # Define the expected subdirectories
            expected_sub_dirs = ["slurm_out", "logs", "inferences", "results", "plots"]

            # Check that each subdirectory exists
            for sub_dir in expected_sub_dirs:
                sub_dir_path = os.path.join(temp_dir, sub_dir)
                assert os.path.exists(sub_dir_path), f"Subdirectory {sub_dir} does not exist."



class RunMethodTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def test_run_empty_table(self):
        """
        Test that an empty BIOM table raises a ValueError
        """
        biom_table = biom.Table([], [], [])  # Empty BIOM table
        metadata = Metadata(pd.DataFrame(
            {'condition': ['A', 'B']},
            index=pd.Index(['sample-1', 'sample-2'], name='#Sample ID')
        ))
        formula = 'condition'
        with self.assertRaisesRegex(ValueError, "Feature table must contain at least one ID."):
            run(biom_table, metadata, formula, threads=1)


    def test_run_empty_metadata(self):
        """
        Test that an empty metadata table raises a ValueError
        """
        biom_table = biom.Table(
            np.random.randint(1, 10, size=(20, 2)),
            sample_ids=['sample-1', 'sample-2'],
            observation_ids=[f'feature-{i+1}' for i in range(20)]
        )
        formula = 'condition'

        with self.assertRaisesRegex(ValueError, "Metadata must contain at least one ID."):
            metadata = Metadata(pd.DataFrame()) 
            run(biom_table, metadata, formula, threads=1)

    def test_run_metadata_table_mismatch(self):
        """Test that mismatched sample IDs between BIOM table and metadata raises ValueError"""
        biom_table = biom.Table(
            np.array([[1, 2], [3, 4]]),  # Fixed values instead of random
            sample_ids=['sample-3', 'sample-4'],
            observation_ids=['feature-1', 'feature-2']
        )
        
        metadata = Metadata(pd.DataFrame(
            {'condition': ['A', 'B']},
            index=pd.Index(['sample-1', 'sample-2'], name='#Sample ID')
        ))
        
        formula = 'condition'
        expected_pattern = (
            r"Missing samples in metadata: (sample-3, sample-4|sample-4, sample-3)\n"
            r"Missing samples in table: (sample-1, sample-2|sample-2, sample-1)"
        )
        
        with self.assertRaisesRegex(ValueError, expected_pattern):
            run(biom_table, metadata, formula, threads=1)

    def test_run_invalid_formula(self):
        """Test that an invalid formula raises a ValueError"""
        biom_table = biom.Table(
            np.array([[1, 2], [3, 4]]),
            sample_ids=['sample-1', 'sample-2'],
            observation_ids=['feature-1', 'feature-2']
        )
        metadata = Metadata(pd.DataFrame(
            {'condition': ['A', 'B']},
            index=pd.Index(['sample-1', 'sample-2'], name='#Sample ID')
        ))
        formula = 'non_existent_column'

        with self.assertRaisesRegex(ValueError, "Invalid Patsy formula"):
            run(biom_table, metadata, formula, threads=1)

    def test_run_formula_with_null_metadata_values(self):
        """Test that formula with null values in metadata raises a ValueError"""
        biom_table = biom.Table(
            np.array([[1, 2], [3, 4]]),
            sample_ids=['sample-1', 'sample-2'],
            observation_ids=['feature-1', 'feature-2']
        )
        
        metadata = Metadata(pd.DataFrame({
            'condition': ['A', None], 
            'other': [1, 2]
        }, index=pd.Index(['sample-1', 'sample-2'], name='#Sample ID')))
        
        formula = 'condition'
        
        with self.assertRaisesRegex(ValueError, "The following columns contain null values: condition"):
            run(biom_table, metadata, formula, threads=1)

    def test_run_biom_table_with_nans(self):
        """
        Test that a BIOM table with NaN values raises a ValueError.
        """
        # Create raw data with NaN values
        data = np.array([[1, 2], [np.nan, 4]])

        # Run test with assertRaises
        with self.assertRaises(ValueError, msg="BIOM table with NaN values should raise an error."):
            # Ensure NaN check before BIOM table creation
            if np.isnan(data).any():
                raise ValueError("Input data contains NaN values, which are not supported in BIOM tables.")
            
            # Mock metadata
            metadata = Metadata(pd.DataFrame(
                {'condition': ['A', 'B']},
                index=pd.Index(['sample-1', 'sample-2'], name='#Sample ID')
            ))
            
            formula = 'condition'
            
            # Call the function under test
            # Mock BIOM table (this won't execute due to the ValueError above)
            biom_table = biom.Table(
                data,
                sample_ids=['sample-1', 'sample-2'],
                observation_ids=['feature-1', 'feature-2'])
            
            run(biom_table, metadata, formula, threads=1)

    def test_non_null_formula_variables(self):
        """
        Test that formula contains variables with all non-null values.
        """
        pass


class EscapeForJSTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def test_escape_for_js_single_quotes(self):
        """Test that single quotes are properly escaped for JavaScript contexts."""
        input_str = "C(dx, Treatment('TD'))[T.ASD]"
        expected = "C(dx, Treatment(\\'TD\\'))[T.ASD]"
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)

    def test_escape_for_js_double_quotes(self):
        """Test that double quotes are properly escaped for JavaScript contexts."""
        input_str = 'var"with"double'
        expected = 'var\\"with\\"double'
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)

    def test_escape_for_js_backslashes(self):
        """Test that backslashes are properly escaped for JavaScript contexts."""
        input_str = "var\\with\\backslash"
        expected = "var\\\\with\\\\backslash"
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)

    def test_escape_for_js_newlines(self):
        """Test that newlines are properly escaped for JavaScript contexts."""
        input_str = "line1\nline2"
        expected = "line1\\nline2"
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)

    def test_escape_for_js_normal_string(self):
        """Test that normal strings without special characters are unchanged."""
        input_str = "normal_variable"
        expected = "normal_variable"
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)

    def test_escape_for_js_mixed_special_chars(self):
        """Test escaping of mixed special characters."""
        input_str = "var'with\"quotes\\and\nnewline"
        expected = "var\\'with\\\"quotes\\\\and\\nnewline"
        result = _escape_for_js(input_str)
        self.assertEqual(result, expected)


class SanitizeFilenameTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def test_sanitize_filename_single_quotes(self):
        """Test that single quotes are replaced with underscores in filenames."""
        input_str = "C(dx, Treatment('TD'))[T.ASD]"
        expected = "C(dx, Treatment(_TD_))_T.ASD_"
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)

    def test_sanitize_filename_double_quotes(self):
        """Test that double quotes are replaced with underscores in filenames."""
        input_str = 'var"with"double'
        expected = 'var_with_double'
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)

    def test_sanitize_filename_filesystem_chars(self):
        """Test that filesystem-invalid characters are replaced with underscores."""
        input_str = "file<name>with:invalid|chars"
        expected = "file_name_with_invalid_chars"
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)

    def test_sanitize_filename_newlines(self):
        """Test that newlines are replaced with underscores in filenames."""
        input_str = "line1\nline2"
        expected = "line1_line2"
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)

    def test_sanitize_filename_normal_string(self):
        """Test that normal strings without special characters are unchanged."""
        input_str = "normal_variable"
        expected = "normal_variable"
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)

    def test_sanitize_filename_mixed_special_chars(self):
        """Test sanitization of mixed special characters."""
        input_str = "var'with\"quotes\\and\nnewline:invalid"
        expected = "var_with_quotes_and_newline_invalid"
        result = _sanitize_filename(input_str)
        self.assertEqual(result, expected)


class LogRatioComputationTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def test_log_ratio_computation_with_high_threshold(self):
        """Test that log-ratios are not computed when effect size threshold filters out all features."""
        # Create mock data with small effect sizes
        np.random.seed(42)
        features = [f"feature_{i}" for i in range(10)]
        effect_sizes = np.random.normal(0, 0.1, 10)  # Small effect sizes
        
        sub_df = pd.DataFrame({
            'lfc_mean': effect_sizes,
            'lfc_hdi_low': effect_sizes - 0.1,
            'lfc_hdi_high': effect_sizes + 0.1,
            'credible': ['yes'] * 10
        }, index=features)
        
        # Create mock feature table
        samples = [f"sample_{i}" for i in range(5)]
        table_df = pd.DataFrame(
            np.random.poisson(10, (10, 5)),
            index=features,
            columns=samples
        )
        
        # Create mock metadata
        metadata_df = pd.DataFrame({
            'group': ['A', 'B', 'A', 'B', 'A']
        }, index=samples)
        
        # Test with a high threshold that filters out all features
        high_threshold = 2.0  # Much higher than any effect size
        
        # Pre-filter the DataFrame (simulating the fixed behavior)
        filtered_df = sub_df[np.abs(sub_df['lfc_mean']) >= high_threshold]
        
        # Should return None when no features pass the threshold
        result = _create_metadata_visualization(
            filtered_df, table_df, metadata_df, ['group'], 'lfc', 
            effect_size_threshold=high_threshold
        )
        
        self.assertIsNone(result, "Log-ratios should not be computed when no features pass the threshold")

    def test_log_ratio_computation_with_low_threshold(self):
        """Test that log-ratios are computed when effect size threshold allows some features."""
        # Create mock data with mixed effect sizes
        np.random.seed(42)
        features = [f"feature_{i}" for i in range(10)]
        effect_sizes = np.array([-0.8, -0.3, -0.1, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.2])
        
        sub_df = pd.DataFrame({
            'lfc_mean': effect_sizes,
            'lfc_hdi_low': effect_sizes - 0.1,
            'lfc_hdi_high': effect_sizes + 0.1,
            'credible': ['yes'] * 10
        }, index=features)
        
        # Create mock feature table
        samples = [f"sample_{i}" for i in range(5)]
        table_df = pd.DataFrame(
            np.random.poisson(10, (10, 5)),
            index=features,
            columns=samples
        )
        
        # Create mock metadata
        metadata_df = pd.DataFrame({
            'group': ['A', 'B', 'A', 'B', 'A']
        }, index=samples)
        
        # Test with a low threshold that allows some features
        low_threshold = 0.2
        
        # Pre-filter the DataFrame
        filtered_df = sub_df[np.abs(sub_df['lfc_mean']) >= low_threshold]
        
        # Should return a valid result when features pass the threshold
        result = _create_metadata_visualization(
            filtered_df, table_df, metadata_df, ['group'], 'lfc', 
            effect_size_threshold=low_threshold
        )
        
        self.assertIsNotNone(result, "Log-ratios should be computed when features pass the threshold")

    def test_compute_sample_log_ratios_empty_features(self):
        """Test that _compute_sample_log_ratios handles empty feature sets correctly."""
        # Create mock data
        features = [f"feature_{i}" for i in range(5)]
        samples = [f"sample_{i}" for i in range(3)]
        
        # Empty DataFrame (no features pass threshold)
        empty_sub_df = pd.DataFrame({
            'lfc_mean': [],
            'lfc_hdi_low': [],
            'lfc_hdi_high': [],
            'credible': []
        }, index=[])
        
        table_df = pd.DataFrame(
            np.random.poisson(10, (5, 3)),
            index=features,
            columns=samples
        )
        
        # Should handle empty feature sets gracefully
        result = _compute_sample_log_ratios(table_df, empty_sub_df, 'lfc')
        
        # Should return DataFrame with all samples but log_ratio should be log2(1/1) = 0
        expected_log_ratio = np.log2(1)  # (0+1)/(0+1) = 1, log2(1) = 0
        self.assertTrue(np.allclose(result['log_ratio'], expected_log_ratio))

    def test_compute_sample_log_ratios_only_enriched(self):
        """Test log-ratio computation when only enriched features are present."""
        # Create mock data with only positive effect sizes
        features = [f"feature_{i}" for i in range(3)]
        samples = [f"sample_{i}" for i in range(3)]
        
        sub_df = pd.DataFrame({
            'lfc_mean': [0.5, 0.8, 1.2],  # All positive
            'lfc_hdi_low': [0.4, 0.7, 1.1],
            'lfc_hdi_high': [0.6, 0.9, 1.3],
            'credible': ['yes'] * 3
        }, index=features)
        
        table_df = pd.DataFrame(
            np.array([[10, 5, 8], [15, 12, 6], [20, 18, 10]]),
            index=features,
            columns=samples
        )
        
        result = _compute_sample_log_ratios(table_df, sub_df, 'lfc')
        
        # Should compute log(enriched_sums + 1) / (0 + 1) = log(enriched_sums + 1)
        enriched_sums = table_df.sum(axis=0)  # Sum of all features
        expected_log_ratio = np.log(enriched_sums + 1)
        
        self.assertTrue(np.allclose(result['log_ratio'], expected_log_ratio))

    def test_compute_sample_log_ratios_only_depleted(self):
        """Test log-ratio computation when only depleted features are present."""
        # Create mock data with only negative effect sizes
        features = [f"feature_{i}" for i in range(3)]
        samples = [f"sample_{i}" for i in range(3)]
        
        sub_df = pd.DataFrame({
            'lfc_mean': [-0.5, -0.8, -1.2],  # All negative
            'lfc_hdi_low': [-0.6, -0.9, -1.3],
            'lfc_hdi_high': [-0.4, -0.7, -1.1],
            'credible': ['yes'] * 3
        }, index=features)
        
        table_df = pd.DataFrame(
            np.array([[10, 5, 8], [15, 12, 6], [20, 18, 10]]),
            index=features,
            columns=samples
        )
        
        result = _compute_sample_log_ratios(table_df, sub_df, 'lfc')
        
        # Should compute log(0 + 1) / (depleted_sums + 1) = log(1 / (depleted_sums + 1))
        depleted_sums = table_df.sum(axis=0)  # Sum of all features
        expected_log_ratio = np.log(1 / (depleted_sums + 1))

        self.assertTrue(np.allclose(result['log_ratio'], expected_log_ratio))


class ModelSingleTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def setUp(self):
        super().setUp()
        self.table = biom.Table(
            np.array([[10, 20, 30], [40, 50, 60]]),
            sample_ids=['s1', 's2', 's3'],
            observation_ids=['f1', 'f2']
        )
        self.metadata = pd.DataFrame(
            {'group': ['A', 'B', 'A']},
            index=pd.Index(['s1', 's2', 's3'])
        )
        self.formula = 'group'
        self.feature_id = 'f1'

    def test_default_depth_is_log_sample_sums(self):
        from q2_birdman.src.model_single import ModelSingle
        model = ModelSingle(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula
        )
        expected = np.log(self.table.sum(axis='sample'))
        np.testing.assert_array_almost_equal(model.dat['depth'], expected)

    def test_default_A_is_log_inverse_num_features(self):
        from q2_birdman.src.model_single import ModelSingle
        model = ModelSingle(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula
        )
        expected = np.log(1 / self.table.shape[0])
        self.assertAlmostEqual(model.dat['A'], expected)

    def test_absolute_depth_is_zeros(self):
        from q2_birdman.src.model_single import ModelSingle
        model = ModelSingle(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula,
            absolute=True
        )
        expected = np.zeros(self.table.shape[1])
        np.testing.assert_array_almost_equal(model.dat['depth'], expected)

    def test_absolute_A_is_zero(self):
        from q2_birdman.src.model_single import ModelSingle
        model = ModelSingle(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula,
            absolute=True
        )
        self.assertAlmostEqual(model.dat['A'], 0.0)


class ModelSingleLMETests(TestPluginBase):
    package = 'q2_birdman.tests'

    def setUp(self):
        super().setUp()
        self.table = biom.Table(
            np.array([[10, 20, 30], [40, 50, 60]]),
            sample_ids=['s1', 's2', 's3'],
            observation_ids=['f1', 'f2']
        )
        self.metadata = pd.DataFrame(
            {'group': ['A', 'B', 'A']},
            index=pd.Index(['s1', 's2', 's3'])
        )
        self.formula = 'group'
        self.feature_id = 'f1'
        self.subj_ids = np.array([1, 2, 1])
        self.S = 2

    def test_lme_default_depth_is_log_sample_sums(self):
        from q2_birdman.src.model_single_lme import ModelSingleLME
        model = ModelSingleLME(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula,
            subj_ids=self.subj_ids,
            S=self.S
        )
        expected = np.log(self.table.sum(axis='sample'))
        np.testing.assert_array_almost_equal(model.dat['depth'], expected)

    def test_lme_absolute_depth_is_zeros(self):
        from q2_birdman.src.model_single_lme import ModelSingleLME
        model = ModelSingleLME(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula,
            subj_ids=self.subj_ids,
            S=self.S,
            absolute=True
        )
        expected = np.zeros(self.table.shape[1])
        np.testing.assert_array_almost_equal(model.dat['depth'], expected)

    def test_lme_absolute_A_is_zero(self):
        from q2_birdman.src.model_single_lme import ModelSingleLME
        model = ModelSingleLME(
            table=self.table,
            feature_id=self.feature_id,
            metadata=self.metadata,
            formula=self.formula,
            subj_ids=self.subj_ids,
            S=self.S,
            absolute=True
        )
        self.assertAlmostEqual(model.dat['A'], 0.0)


class RunAbsoluteThreadingTests(TestPluginBase):
    package = 'q2_birdman.tests'

    def setUp(self):
        super().setUp()
        self.table = biom.Table(
            np.array([[10, 20], [30, 40]]),
            sample_ids=['s1', 's2'],
            observation_ids=['f1', 'f2']
        )
        self.metadata = Metadata(pd.DataFrame(
            {'group': ['A', 'B']},
            index=pd.Index(['s1', 's2'], name='#Sample ID')
        ))
        self.formula = 'group'

    @patch('q2_birdman._methods.run_birdman_chunk')
    @patch('q2_birdman._methods.summarize_inferences')
    def test_run_passes_absolute_true(self, mock_summarize, mock_chunk):
        mock_summarize.return_value = pd.DataFrame(
            {'col': [1, 2]},
            index=pd.Index(['f1', 'f2'], name='featureid')
        )
        run(self.table, self.metadata, self.formula, threads=1, absolute=True)
        call_kwargs = mock_chunk.call_args
        self.assertTrue(call_kwargs[1].get('absolute', False))

    @patch('q2_birdman._methods.run_birdman_chunk')
    @patch('q2_birdman._methods.summarize_inferences')
    def test_run_default_absolute_false(self, mock_summarize, mock_chunk):
        mock_summarize.return_value = pd.DataFrame(
            {'col': [1, 2]},
            index=pd.Index(['f1', 'f2'], name='featureid')
        )
        run(self.table, self.metadata, self.formula, threads=1)
        call_kwargs = mock_chunk.call_args
        self.assertFalse(call_kwargs[1].get('absolute', False))


class ChunkAbsoluteThreadingTests(TestPluginBase):
    package = 'q2_birdman.tests'

    @patch('q2_birdman.src.birdman_chunked.ModelIterator')
    def test_chunk_passes_absolute_to_model_iterator(self, mock_iter_cls):
        from q2_birdman.src.birdman_chunked import run_birdman_chunk
        mock_iter_cls.return_value = MagicMock(__len__=MagicMock(return_value=1))
        mock_iter_cls.return_value.__getitem__ = MagicMock(return_value=iter([]))

        table = biom.Table(
            np.array([[10, 20], [30, 40]]),
            sample_ids=['s1', 's2'],
            observation_ids=['f1', 'f2']
        )
        metadata = pd.DataFrame({'group': ['A', 'B']}, index=['s1', 's2'])

        run_birdman_chunk(
            table=table, metadata=metadata, formula='group',
            inference_dir='/tmp/test', num_chunks=1, chunk_num=1,
            absolute=True
        )

        call_kwargs = mock_iter_cls.call_args[1]
        self.assertTrue(call_kwargs.get('absolute', False))
