from utils import *

__all__ = ['RNASuiteDeseqParameters', 'RNASuiteShowDEGParameters',
	'RNASuiteEdgerParameters', 'RNASuiteDEGDistPlotParameters',
	'RNASuiteDEGVolcanoPlotParameters', 'RNASuiteDEGVennPlotParameters',
	'RNASuiteDEGUpsetPlotParameters', 'RNASuiteDeseqMaPlotControlParameters'
]

class RNASuiteParameter(dict):
	def __init__(self, expose=True, default=None, index=False, options=[], **kwargs):
		self.update(kwargs)
		self['expose'] = expose
		self['default'] = default

		if kwargs['type'] == 'list':
			self['index'] = index
			self['options'] = options

		elif kwargs['type'] == 'select':
			self['options'] = options

	def __getattr__(self, key):
		return self[key]

class RNASuiteParameters:
	def __init__(self, *args):
		self.params = args
		self.mapping = {p.key: i for i, p in enumerate(self.params)}

	def __getitem__(self, key):
		return self.params[self.mapping[key]]

	def __contains__(self, key):
		return key in self.mapping

	def __iter__(self):
		for p in self.params:
			yield p

RNASuiteDeseqParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'fdr',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		decimal = 5,
		display = 'FDR:',
		default = 0.05
	),
	RNASuiteParameter(
		key = 'logfc',
		type = 'int',
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	RNASuiteParameter(
		key = 'compare',
		type = 'list',
		display = 'Comparison between:'
	),
	RNASuiteParameter(
		key = 'control',
		type = 'list',
		display = 'Control group:'
	),
	RNASuiteParameter(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
	),
	RNASuiteParameter(
		key = 'eliminate',
		type = 'select',
		display = 'Consider effect of factors:',
		expose = False
	),
	RNASuiteParameter(
		key = 'design',
		type = 'str',
		display = "Model design formula:"
	),
	RNASuiteParameter(
		key = 'custom',
		type = 'bool',
		display = "Custom design formula:",
		expose = False
	)
)

RNASuiteEdgerParameters = RNASuiteParameters(
	AttrDict(
		key = 'fdr',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = 'int',
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'compare',
		type = 'list',
		options = [],
		display = 'Comparison between:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'control',
		type = 'list',
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'eliminate',
		type = 'select',
		options = [],
		display = 'Considering effect of factors:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'design',
		type = 'str',
		display = "Model design formula:",
		default = None
	),
	AttrDict(
		key = 'custom',
		type = 'bool',
		display = "Custom design formula:",
		default = False
	),
	AttrDict(
		key = 'replicate',
		type = 'list',
		options = ['Biological Replicates', 'No Biological Replicates'],
		display = "Use edgeR with:",
		default = None,
		index = True
	),
	AttrDict(
		key = 'bcv',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		display = "BCV value:",
		default = 0.4
	),
	AttrDict(
		key = 'method',
		type = 'list',
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
		type = 'float',
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = 'int',
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'control',
		type = 'list',
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	)
)

RNASuiteDEGDistPlotParameters = RNASuiteParameters(
	AttrDict(
		key = 'contrasts',
		type = 'contrast',
		default = None,
		display = 'Contrasts:'
	),
	AttrDict(
		key = 'plot',
		type = 'list',
		options = [
			'Stacked bar plot',
			'Dodged bar plot'
		],
		display = "Plot type:",
		default = None,
		index = True
	),
	AttrDict(
		key = 'showval',
		type = 'bool',
		display = "Show value labels:",
		default = False
	),
	AttrDict(
		key = 'rotate',
		type = 'int',
		range = [0, 90],
		step = 5,
		display = "X labels rotate angle:",
		default = 0
	),
	AttrDict(
		key = 'ucolor',
		type = 'color',
		display = "Up-regulated color:",
		default = '#e41a1c'
	),
	AttrDict(
		key = 'dcolor',
		type = 'color',
		display = "Down-regulated color:",
		default = '#377eb8'
	)
)

RNASuiteDEGVolcanoPlotParameters = RNASuiteParameters(
	AttrDict(
		key = 'ncolor',
		type = 'color',
		display = "Not-significant color:",
		default = '#eaeded'
	),
	AttrDict(
		key = 'ucolor',
		type = 'color',
		display = "Up-regulated color:",
		default = '#e41a1c'
	),
	AttrDict(
		key = 'dcolor',
		type = 'color',
		display = "Down-regulated color:",
		default = '#377eb8'
	),
	AttrDict(
		key = 'line',
		type = 'bool',
		display = "Add threshold dashed lines:",
		default = True
	),
	AttrDict(
		key = 'top',
		type = 'int',
		range = (0, 1000),
		step = 1,
		display = "Show top significant gene names:",
		default = 10
	),
	AttrDict(
		key = 'gname',
		type = 'list',
		options = ('gene ID column', 'gene annotation'),
		display = 'Show gene names from:',
		default = None,
		index = True
	),
	AttrDict(
		key = 'gnsep',
		type = 'str',
		display = "Gene ID separator:",
		default = '|'
	),
	AttrDict(
		key = 'gncol',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 1,
		display = "Show gene names in column:"
	),
	AttrDict(
		key = 'limit',
		type = 'int',
		range = (0, 100),
		step = 1,
		default = 0,
		display = "Log2 Fold Change limit:"
	)
)

RNASuiteDEGVennPlotParameters = RNASuiteParameters(
	AttrDict(
		key = 'contrasts',
		type = 'contrast',
		default = None,
		display = 'Contrasts:'
	),
	AttrDict(
		key = 'colors',
		type = 'colors',
		default = None,
		display = 'Contrast colors:'
	),
	AttrDict(
		key = 'opacity',
		type = 'float',
		range = (0, 1),
		step = 0.1,
		default = 0.5,
		display = "Color opacity:"
	),
	AttrDict(
		key = 'percent',
		type = 'bool',
		display = "Show percent values:",
		default = False
	),
	AttrDict(
		key = 'degtype',
		type = 'list',
		display = "DEG Type:",
		options = ('All DEGs', 'Up-regulated', 'Down-regulated'),
		default = None,
		index = True
	)
)

RNASuiteDEGUpsetPlotParameters = RNASuiteParameters(
	AttrDict(
		key = 'contrasts',
		type = 'contrast',
		default = None,
		display = 'Contrasts:'
	),
	AttrDict(
		key = 'degtype',
		type = 'list',
		display = "DEG Type:",
		options = ('All DEGs', 'Up-regulated', 'Down-regulated'),
		default = None,
		index = True
	)
)

RNASuiteDeseqMaPlotControlParameters = RNASuiteParameters(
	AttrDict(
		key = 'main',
		type = 'str',
		default = '',
		display = "Title:"
	),
	AttrDict(
		key = 'xlab',
		type = 'str',
		default = 'mean of normalized counts',
		display = "X label:"
	),
	AttrDict(
		key = 'ylim',
		type = 'int',
		range = (-1000, 1000),
		step = 1,
		default = 0,
		display = "X Limits:"
	),
	AttrDict(
		key = 'colNonSig',
		type = 'color',
		default = 'gray60',
		display = "Non-significant point color:"
	),
	AttrDict(
		key = 'colSig',
		type = 'color',
		default = 'blue',
		display = "Significant point color:"
	),
	AttrDict(
		key = 'colLine',
		type = 'color',
		default = 'gray40',
		display = "Horizontal line color:"
	),
)
