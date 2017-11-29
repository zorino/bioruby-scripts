#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-11-01
# version: 	0.01

require 'bio'

def get_species species_file

  species = []
  File.open(species_file) do |f|
    while l = f.gets
      s, nb_of_genomes = l.chomp!.split("\t")
      species.push([s, nb_of_genomes.to_i])
    end
  end
  species_sorted = species.sort_by{ |k,v| v }.reverse

  return species_sorted

end


def extract_consensus_seq rna, matrix_f, fasta_f

  canditate = []
  best_rna = rna.sort_by { |x,y,z| y }.reverse[0]
  index = best_rna[0].gsub("cluster_","").to_i
  iterator = -1
  File.open(matrix_f, "r") do |f|
    f.gets                      # skip header
    while l = f.gets
      iterator += 1
      break if iterator > index
      next if iterator != index
      # in THE cluster
      lA = l.chomp!.split("\t")
      lA.each do |x|
        canditate.push(*x.split(";"))
      end
    end
  end

  consensus = Bio::Alignment.new()
  flat_f = Bio::FlatFile.auto(fasta_f)
  flat_f.each_entry do |e|
    if canditate.include?(e.definition) and e.seq.length == best_rna[2]
      consensus.add_seq(e.seq, e.definition)
    end
  end

  consensus_seq = []
  (0..best_rna[2]-1).each do |i|
    pos = []
    consensus.each_seq do |s|
      pos << s[i]
    end
    freq = pos.inject(Hash.new(0)) { |h,v| h[v] += 1; h }
    pos_nt = pos.max_by { |v| freq[v] }
    consensus_seq << pos_nt
  end

  return consensus_seq.join('')

end

# Will extract a consensus for all RNA sequences with the same product name
# used for RNAs but not for CDS since different protein can have the same name (hypothetical protein)
def extract_rna_consensus_cdhit cdhit_f, species_name

  fasta_f = cdhit_f.gsub(".cd-hit", "")
  clstr_f = cdhit_f + ".clstr"
  info_f = clstr_f + ".info.tsv"
  matrix_f = clstr_f + ".matrix.tsv"
  path = File.dirname(fasta_f)

  rnas = {}
  File.open(info_f, "r") do |f|
    f.gets
    while l = f.gets
      cluster_id, occurences, length, name, others = l.chomp!.split("\t")
      if ! rnas.has_key? name
        rnas[name] = []
      end
      rnas[name].push([cluster_id, occurences.split("/")[0].to_i, length.to_i])
    end
  end

  if ! Dir.exists? "#{path}/zz_consensus"
    Dir.mkdir "#{path}/zz_consensus"
  end

  rnas.each do |k, v|
    consensus_seq = extract_consensus_seq v, matrix_f, fasta_f
    if consensus_seq != ""
      bioseq = Bio::Sequence.new(consensus_seq)
      # puts bioseq.output_fasta("#{species_name} #{k}")
      output_file = k.gsub(" ","_") + ".fasta"
      File.open("#{path}/zz_consensus/#{output_file}", "w") do |f|
        f.write(bioseq.output_fasta("#{species_name} #{k}"))
      end
    end
  end

  return ""

end


# Extract CDS can take several CDS with same product name
def extract_cds_consensus_cdhit cdhit_f, species_name

  fasta_f = cdhit_f.gsub(".cd-hit", "")
  clstr_f = cdhit_f + ".clstr"
  info_f = clstr_f + ".info.tsv"
  matrix_f = clstr_f + ".matrix.tsv"
  path = File.dirname(fasta_f)

  # cluster info
  clusters = {}
  genome_origins = []
  File.open(info_f, "r") do |f|
    f.gets
    while l = f.gets
      cluster_id, occurences, length, name, others = l.chomp!.split("\t")
      if name == "chromosomal replication initiation protein DnaA"
        genome_origins.push(cluster_id)
      end
      if ! clusters.has_key? name
        clusters[cluster_id] = []
      end
      clusters[cluster_id].push([name, occurences.split("/")[0].to_i, length.to_i])
    end
  end


  # genomes matrix
  genomes = {}
  genomes_list = []
  File.open(matrix_f, "r") do |f|

    # get header genomes
    l = f.gets
    l.chomp!.split("\t")[1..-1].each do |g|
      genomes_list.push(g)
      genomes[g] = {}
      genomes[g]['all_cds'] = []
      genomes[g]['origin_cluster'] = ""
    end

    p genome_origins

    while l = f.gets
      lA = l.chomp!.split("\t")
      cluster_id = lA[0]
      lA[1..-1].each_with_index do |g,i|
        last_origin = ""
        g.split(";").each do |hit|
          next if hit == "-"
          sample, genome_id, prot_id, prot_locus, gene, product = hit.split("|")
          genomes[genomes_list[i]][prot_locus + "_-_" + genome_id] = cluster_id
          genomes[genomes_list[i]]['all_cds'].push(prot_locus + "_-_" + genome_id)
          if genome_origins.include? cluster_id
            if last_origin == ""
              genomes[genomes_list[i]]['origin_cluster'] = cluster_id
              genomes[genomes_list[i]]['origin_cds'] = prot_locus + "_-_" + genome_id
              last_origin = prot_locus
            elsif last_origin > prot_locus
              genomes[genomes_list[i]]['origin_cluster'] = cluster_id
              genomes[genomes_list[i]]['origin_cds'] = prot_locus + "_-_" + genome_id
            end
          end
        end

      end

    end

  end

  # sort genomes cds
  genomes.each_key do |g|
    genomes[g]['all_cds'].sort!
    # first = genomes[g]['all_cds'][0]
    # puts first + " " + genomes[g][first]
  end

  puts "\n# Results :"
  # get genome origin
  genomes.each_key do |g|
    puts genomes[g]['origin_cluster']
    puts genomes[g]['origin_cds']
  end


  # walking for consensus order from clusters


end



usage = "

create-species-model-sum.rb <species name> <species dir> <seq type>

  seq_type = rRNA, tRNA, CDS

"

if ARGV.length < 3
  abort usage
end

species_name = ARGV[0]
species_dir = ARGV[1]
seq_type = ARGV[2]

# species = get_species species_file

if seq_type.downcase == "rrna"

  rna_cdhit_file = species_dir+"/zz_pangenome/rRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name

elsif seq_type.downcase == "trna"

  rna_cdhit_file = species_dir+"/zz_pangenome/tRNA.fasta.cd-hit"
  extract_rna_consensus_cdhit rna_cdhit_file, species_name

elsif seq_type.downcase == "cds"

  cds_cdhit_file = species_dir+"/zz_pangenome/CDS.fasta.cd-hit"
  extract_cds_consensus_cdhit cds_cdhit_file, species_name

end

