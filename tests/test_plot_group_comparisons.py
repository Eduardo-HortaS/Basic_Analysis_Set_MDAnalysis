import os
import pickle
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plotting import plot_group_comparisons as pgc


def _write_rmsf_pickle(path, selection='protein and backbone', ref_selection='chainid A'):
    data = {
        'rmsf': [1.0, 2.0],
        'resids': [1, 2],
        'selection': selection,
        'ref_selection': ref_selection,
    }
    with open(path, 'wb') as f:
        pickle.dump(data, f)


class TestPlotRmsfGroups:
    def test_uses_pickle_metadata_for_selection_label(self, tmp_path, monkeypatch):
        work_dir = tmp_path / 'work'
        work_dir.mkdir()
        output_dir = tmp_path / 'plots'
        output_dir.mkdir()

        _write_rmsf_pickle(work_dir / 'rmsf_plot_SysA_v1_rep1.pkl')
        _write_rmsf_pickle(work_dir / 'rmsf_plot_SysB_v1_rep1.pkl')

        captured = {}

        def fake_plot_rmsf_comparison_average(pickle_groups, labels, gname, **kwargs):
            captured['gname'] = gname
            captured['labels'] = labels
            captured['selection_label'] = kwargs.get('selection_label')

        monkeypatch.setattr(pgc, '_import_rmsf_plotters', lambda: (None, fake_plot_rmsf_comparison_average))

        produced = pgc._plot_rmsf_groups(
            {'grp': [('SysA', 'v1'), ('SysB', 'v1')]},
            work_dir=str(work_dir),
            output_dir=str(output_dir),
            num_replicates=1,
            dpi=200,
            replicate_mode='average',
            target_selection='fallback target',
            ref_selection='fallback ref',
        )

        assert produced == 1
        assert captured['gname'] == 'grp'
        assert captured['labels'] == ['SysA / v1', 'SysB / v1']
        assert captured['selection_label'] == 'target=protein and backbone | ref=chainid A'
