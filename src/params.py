from utils import *

__all__ = ['RNASuiteDEGParameters', 'RNASuiteShowDEGParameters',
	'RNASuiteEdgerParameters'
]

class RNASuiteParameters():
	def __init__(self, *args):
		self.params = args
		self.mapping = {p.key: i for i, p in enumerate(self.params)}

	def __getitem__(self, key):
		return self.params[self.mapping[key]]

	def __iter__(self):
		for p in self.params:
			yield p

RNASuiteDEGParameters = RNASuiteParameters(
	AttrDict(
		key = 'fdr',
		type = float,
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = int,
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'compare',
		type = list,
		options = [],
		display = 'Comparison between:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'control',
		type = list,
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = list,
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'eliminate',
		type = set,
		options = [],
		display = 'Eliminate effect of factors:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'design',
		type = str,
		display = "Model design formula:",
		default = None
	),
	AttrDict(
		key = 'custom',
		type = bool,
		display = "Custom design formula:",
		default = False
	)
)

RNASuiteEdgerParameters = RNASuiteParameters(
	AttrDict(
		key = 'fdr',
		type = float,
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = int,
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'compare',
		type = list,
		options = [],
		display = 'Comparison between:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'control',
		type = list,
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = list,
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'eliminate',
		type = set,
		options = [],
		display = 'Considering effect of factors:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'design',
		type = str,
		display = "Model design formula:",
		default = None
	),
	AttrDict(
		key = 'custom',
		type = bool,
		display = "Custom design formula:",
		default = False
	),
	AttrDict(
		key = 'replicate',
		type = list,
		options = ['Biological Replicates', 'No Biological Replicates'],
		display = "Use edgeR with:",
		default = None,
		index = True
	),
	AttrDict(
		key = 'bcv',
		type = float,
		range = (0, 1),
		step = 0.01,
		display = "BCV value:",
		default = 0.4
	),
	AttrDict(
		key = 'method',
		type = list,
		options = [
			"Quasi-likelihood F-tests",
			"Likelihood ratio tests",
			"Classic exact test"
		],
		display = "Testing method:",
		default = None,
		index = True
	)
)

RNASuiteShowDEGParameters = RNASuiteParameters(
	AttrDict(
		key = 'fdr',
		type = float,
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = int,
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'control',
		type = list,
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = list,
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	)
)
