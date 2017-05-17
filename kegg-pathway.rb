#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-11-27
# version: 	0.01
# licence:  	

usage = "kegg-pathway.rb <organism>"

require 'bio'

# org = ARGV[0] or abort usage

serv = Bio::KEGG::API.new

# list = serv.list_pathways("eco")
# list.each do |path|
#   print path.entry_id, "\t", path.definition, "\n"
# end


# org= serv.list_organisms
# org.each do |o|
#   print o.entry_id, "\t", o.definition, "\n"
# end

# list=serv.list_pathways(org)
# list.each do |path|
#   print path.entry_id, "\t", path.definition, "\n"
#   genes=serv.get_genes_by_pathway(path.entry_id)
#   genes.each do |g|
#     puts g
#   end
# end
