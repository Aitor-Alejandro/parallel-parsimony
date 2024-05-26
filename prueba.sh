#!/bin/bash
for i in {1..20}; do
	./PARS ../Datasets/neotrop/reference.phy ../Datasets/neotrop/query100.phy ../Datasets/neotrop/reftree.tree g
done
