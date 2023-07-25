import copy
from io import StringIO

from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PhyloXML import Phylogeny

#%matplotlib inline

tree = Phylo.read ("/home/aitor-pc/Escritorio/[TFG]PhyloPlacement/Datasets/salmonella/reftree.tree", "newick")
#with io.open('/home/aitor-pc/Escritorio/[TFG]PhyloPlacement/Datasets/salmonella/reftree.tree', encoding='utf8') as fp:
#    tree = load(fp)
Phylo.draw(tree)
print("Hello world")
