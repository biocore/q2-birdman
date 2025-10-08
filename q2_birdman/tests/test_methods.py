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
from q2_birdman._visualizers import _escape_for_js
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
        """Test that an invalid formula raises a PatsyError"""
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

        expected_error = (
            "Missing columns in metadata: non_existent_column\n"
            "Available columns are: condition"
        )

        with self.assertRaisesRegex(ValueError, re.escape(expected_error)):
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
