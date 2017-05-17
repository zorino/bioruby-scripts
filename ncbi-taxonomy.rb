#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2015-06-04
# version: 	0.01

require 'bio'



ncbi = Bio::NCBI::REST.new
Bio::NCBI.default_email = "default@default.com"

# Acetothermia_bacterium_SCGC_AAA255-C06
# Acholeplasma_axanthum
# Acholeplasma_brassicae
# Acholeplasma_equifetale
# Acholeplasma_granularum
# Acholeplasma_hippikon
# Acholeplasma_laidlawii
# Acholeplasma_modicum
# Acholeplasma_multilocale
# Acholeplasma_oculi
# Acholeplasma_palmae
# Achromobacter_arsenitoxydans
# Achromobacter_insuavis
# Achromobacter_piechaudii
# Achromobacter_sp._DH1f
# Achromobacter_sp._LC458
# Achromobacter_sp._RTa
# Achromobacter_xylosoxidans
# Acidaminococcus_fermentans
# Acidaminococcus_intestini
# Acidaminococcus_sp._BV3L6
# Acidaminococcus_sp._D21
# Acidaminococcus_sp._HPA0509
# Acidicaldus_organivorans



res = ncbi.esearch("Acidaminococcus fermentans",{"db"=>"taxonomy"})
puts res
taxon = ncbi.efetch(res,{"db"=>"taxonomy", "rettype"=>"xml"})
puts taxon

# tax = Bio::Taxonomy.new(taxon)
# puts tax.rank
