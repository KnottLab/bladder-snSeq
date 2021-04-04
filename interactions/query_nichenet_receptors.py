#!/usr/bin/env python

import numpy as np
import pandas as pd
import copy
import os

import argparse
from utils import make_logger



def build_receptor_regulator_lists(receptors, 
                                   weighted_lr_sig, 
                                   max_n=10,
                                   tfs_list=None):
  receptor_regulators_weighted = {}
  for r in sorted(receptors):
    if r.startswith('OR') or r.startswith('TAS'):# or r.startswith('RX'):
      logger.warning(f'Receptor hit a filtered pattern: {r}')
      continue
        
    if r not in weighted_lr_sig['from'].values:
      logger.warning(f'Receptor not in weighted signaling graph: {r}')
      continue
        
    regulators = weighted_lr_sig.loc[weighted_lr_sig['from'] == r, :]
    w = np.array(regulators['weight'])
        
    names = np.array(regulators['to'])
    try: # Clobber names to the shape we need, if it's just one thing.
      n_names = len(names)
    except:
      names = np.expand_dims(names, -1)
    
  #     # Filter things that are also transcription factors
  #     is_tf = np.array([name in transcription_factors for name in names])
  #     if np.sum(is_tf) == 0:
  #         print(f'WARN: {r} no annotated TF partners')
  #     elif np.sum(is_tf) == len(names):
  #         pass
  #     else:
  #         names = names[is_tf]
  #         w = w[is_tf]
  #         print(f'INFO: {r} associated TFs', len(names))
        
    # print(names)
    if len(names) == 1:
      pass
    else:
      order = np.argsort(-w)
      names = names[order]
      w = w[order]
      names = names[:min(len(names), max_n)]

      ## Apply a weight cutoff to subset long lists -- this is bad?
  #         if len(names) > max_n:
  #             w_hi = np.max(w)
  #             w_cutoff = w_hi * 0.5
  #             # Restrict to 50% of top weight
  #             names = names[w > w_cutoff]
  #             w = w[w > w_cutoff]
  #             # enforce the number cutoff
  #             names = names[:max_n]
  #         else:
  #             w = 0
  #     except:
  #         names = np.expand_dims(names, -1).tolist()
  #         print(f'{r}: {names}, {len(names)} {type(names)}')
    
    u_names = np.unique(names)
    logger.info(f'{r} assigning {len(u_names)} signaling partners')
    receptor_regulators_weighted[r] = u_names

  return receptor_regulators_weighted




def query_network(gene, network, top_n_per_step=10):
  rows = network.loc[network['from'] == gene, :]
  w = rows.loc[:, 'weight'].values
  srt = np.argsort(w)[::-1]
  genes = genes[srt][:top_n_per_step]
  return genes


def get_downstream_genes(starts, 
                         network, 
                         steps=1, 
                         edge_threshold=None, 
                         edge_percentile=None, 
                         top_n_per_step=None,
                         top_n_decay=lambda s,t: max(1, t-s)
                        ):
  """ Given a list of starting genes query a network for downstream nodes
  :param starts: a list
  """
  if isinstance(starts, list):
    starts = starts
  else:
    starts = [starts]

  input_starts = starts # For later, because we modify `starts` in place.
  downstream_genes = starts # Downstream genes includes the starting genes
  must_threshold = (edge_threshold is not None or \
                    edge_percentile is not None or \
                    top_n_per_step is not None)
  
  logger.info(f'Searching GRN for downstream genes from {starts}')

  # # Traverse the graph a number of times
  # for sg in input_starts:

  #   if sg not in network['from']:
  #     logger.warning(f'Start gene {sg} not in network from column')
  #     continue

  #   step_targets = query_network(sg, network, top_n_per_step=top_n_per_step)
  #   logger.info(f'gene {sg} targets {step_targets}')
  #   for tg in step_targets:
  #     if tg not in downstream_genes:
  #       downstream_genes.append(tg)

  # return downstream_genes

  for step in range(steps):
    # At each step have a list of starting points and add the first-degree connections to the list
    logger.info(f'step {step}')
    step_genes = []
    for g in starts:
      skip_thresholding = False
      # For each start find entries in the network with this gene as the emitting node
      if g not in network['from'].values:
        logger.warning(f'Start gene {g} not in network from column')
        continue
          
      rows = network.loc[network['from'].values == g, :]
      genes = rows['to']
      
      try: 
        genes = genes.values
      except:
        genes = [genes]
        # If we only get one result, we're not going to toss it 
        skip_thresholding = True
      

      # # Optionally impose conditions on the targets based on the edge weights
      # if rows.shape[0] > 0 and must_threshold and not skip_thresholding:
        
      #   if edge_percentile:
      #     pct = np.quantile(w, edge_percentile)
      #     genes = genes[w >= pct]
      #     w = w[w >= pct]
        
      #   if edge_threshold:
      #     genes = genes[w > edge_threshold]
      #     w = w[w > edge_threshold]
        
      if top_n_per_step:
        w = rows['weight'].values
        srt = np.argsort(w)[::-1]
        genes = genes[srt][:top_n_decay(step, top_n_per_step)]

      logger.info(f'Receptor {step} {g} :--> {genes}')

      step_genes += list(genes)
        
        
    # Avoid tracking a multiplicity of the same gene (we unique them at the end anyway)
    step_genes = list(np.unique(step_genes))
    
    # Avoid repetitively traversing cycles
    step_genes = [g for g in step_genes if g not in downstream_genes]
    
    # Record the genes we picked up this round
    downstream_genes += step_genes

    # Set a new set of starting points for the next iteration
    starts = step_genes
          
  if len(downstream_genes) == 1:
    logger.warning('Downstream genelist from {input_starts} length 1')
  
  downstream_genes = list(np.unique(downstream_genes))
  return downstream_genes
    

  

def get_genes(receptors, receptor_regulators_weighted, weighted_gr, steps=2, top_n_per_step=10):
  receptor_genelists = {}
  for r in sorted(receptors):
    if r not in receptor_regulators_weighted.keys():
      continue
        
    start_genes = copy.copy(receptor_regulators_weighted[r])
    try:
      # Try to make it a list
      start_genes = list(start_genes)
    except:
      # Fails if start_genes is just one string
      start_genes = start_genes

    logger.info(f'Receptor {r} starting GRN traversal with {start_genes}')
    downstream_genes = get_downstream_genes(start_genes, 
                            weighted_gr, 
                            steps=steps, 
                            edge_threshold=None, 
                            edge_percentile=None, 
                            top_n_per_step=top_n_per_step, 
                            top_n_decay=lambda s,t: max(1, t - 2*s))
    # Add the receptor to the list
    downstream_genes += [r]
    receptor_genelists[r] = downstream_genes

  return receptor_genelists




def main(ARGS):

  weighted_lr_sig = pd.read_csv(ARGS.weighted_lr_sig, index_col=0, header=0)
  weighted_gr = pd.read_csv(ARGS.weighted_gr, index_col=0, header=0)
  logger.info(f'Original weighted lr sig table: {weighted_lr_sig.shape}')

  # Intersect signaling partners that are also in the from col of the GRN
  weighted_lr_sig = weighted_lr_sig.loc[weighted_lr_sig['to'].isin(weighted_gr['from'].unique()), :]

  logger.info(f'Trimmed weighted lr sig table: {weighted_lr_sig.shape}')
  logger.info(f'Weighted gr table: {weighted_gr.shape}')

  receptors = [l.strip() for l in open(ARGS.receptors_fname)]
  logger.info(f'Loaded list of {len(receptors)} receptors')

  if not os.path.isdir(ARGS.output_dir):
    logger.info(f'Creating output directory {ARGS.output_dir}')
    os.makedirs(ARGS.output_dir)

  receptor_regulators_weighted = build_receptor_regulator_lists(receptors, 
                                   weighted_lr_sig, 
                                   max_n=ARGS.n_sig_partners,
                                   tfs_list=None)
  logger.info(f'Receptor - regulator lists: {len(receptor_regulators_weighted)}')
  receptor_genelists = get_genes(receptors, receptor_regulators_weighted, weighted_gr,
                                 steps=ARGS.steps)
  
  logger.info(f'Continuing with {len(receptor_genelists)} gene lists')
  
  n_output = 0
  for r, v in sorted(receptor_genelists.items()):
    if len(v) < 10:
      logger.warning(f'Receptor {r} gene list failed length cutoff ({len(v)})')
      continue
    with open(f'{ARGS.output_dir}/{r}.txt', 'w+') as f:
      for v_ in v:
        f.write(f'{v_}\n')
    n_output += 1

  logger.info(f'Recorded {n_output} gene lists from {len(receptors)} input receptors')


if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument(
    'receptors_fname', default=None, type=str,
    help='newline-delimited list of receptors to search'
  )
  parser.add_argument(
    'output_dir', default=None, type=str,
    help='someplace to stash the results. it will be created in case it does not exist.'
  )
  parser.add_argument(
    '--weighted_lr_sig', default='../data/nichenet_weighted_lr_sig.csv', type=str,
    help='a three-column table with from, to, and weight columns representing directed connections '+\
         'between proteins'
  )
  parser.add_argument(
    '--weighted_gr', default='../data/nichenet_weighted_gr.csv', type=str,
    help='a three-column table with from, to, and weight columns representing gene regulatory relationships'
  )
  parser.add_argument(
    '--steps', default=2, type=int,
    help='the number of moves away from each start node to consider. default=2'
  )
  parser.add_argument(
    '--n_sig_partners', default=10, type=int,
    help='the maximum number of elements to take at each step. default=10'
  )

  ARGS = parser.parse_args()

  logger = make_logger()
  for k, v in ARGS.__dict__.items():
    logger.info(f'Argument: {k}: {v}')

  main(ARGS)