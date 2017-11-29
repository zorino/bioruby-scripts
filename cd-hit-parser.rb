#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:
# date:    	14-02-23
# version: 	0.01
# licence:

require 'bio'


class CdhitParser

  def initialize cdHitFile, fastaFile
    @fastaFile = fastaFile
    @cdHitFile = cdHitFile
  end


  def build_genome_index

    # Get Genome Index
    puts "Building Genome Index.."
    genome_index = []
    genome_genes = {}
    fasta_flat = Bio::FlatFile.auto(@fastaFile)
    fasta_flat.each_entry do |f|
      f_info = f.definition.split("|")
      genome_index << f_info[0] if ! genome_index.include? f_info[0]
      genome_genes[f_info[0]] = {} if !genome_genes.has_key? f_info[0]
      genome_genes[f_info[0]][f.definition.split(" ")[0]] = f.definition
    end
    return genome_index, genome_genes

  end


  def parse_cdhit_file

    clusters = []
    cluster_array = {
      id: 0,
      hits: [],
      lengths: [],
      gene_names: []
    }
    cluster_id = -1

    puts "Parsing CD-Hit file.."
    File.open(@cdHitFile,"r") do |f|
      while l=f.gets
        if l[0] == ">"              # cluster id / new cluster
          #   if lineArray[0] != "-"
          #     puts lineArray.join("\t")
          #   end
          if cluster_array[:hits] != []
            clusters << cluster_array.dup
          end
          cluster_id =  l.chomp![1..-1].split(" ")[1].to_i
          cluster_array = nil
          cluster_array = {
            id: cluster_id,
            hits: [],
            lengths: [],
            gene_names: []
          }
        else                        # member of cluster
          lA = l.chomp!.split(", >")
          cluster_array[:lengths].push(lA[0].split()[1].gsub("aa,","").to_i)
          hit = lA[1].split("... ")
          idPercentage = ""
          if ! hit[1].match(/\*/).nil?
            idPercentage = "100"
            # cluster_array.unshift("#{hit[0]}|#{length}|*")
            # cluster_array[:hits].push("#{hit[0].split('|')[0]}")
            cluster_array[:hits].push("#{hit[0]}")
          else
            idPercentage = hit[1].gsub("at ","").gsub("%","")
            # cluster_array.push("#{hit[0]}|#{length}|#{idPercentage}")
            # cluster_array[:hits].push("#{hit[0].split('|')[0]}")
            cluster_array[:hits].push("#{hit[0]}")
          end
        end
      end
      clusters << cluster_array.dup
    end

    return clusters

  end


  def build_cluster_matrix clusters, genome_index, genome_genes

    puts "Building Cluster Matrix.."
    # Clusters-Matrix.tsv
    file_out = File.open(@cdHitFile+".matrix.tsv","w")
    file_out.write("Clusters\t" + genome_index.join("\t") + "\n")
    i = 1
    clusters.each do |c|
      cluster = Array.new(genome_index.size, "-")
      c[:hits].each do |s|
        genome = s.split("|")[0]
        gene = genome_genes[genome][s]
        if cluster[genome_index.index(genome)] != "-"
          cluster[genome_index.index(genome)] = cluster[genome_index.index(genome)] + ";" + gene
        else
          cluster[genome_index.index(genome)] = gene
        end
        c[:gene_names].push(gene.split("|")[-1])
      end
      c[:occurences] = cluster.size - cluster.count("-")
      file_out.write("cluster_#{c[:id]}\t" + cluster.join("\t") + "\n")
      cluster = nil
      i += 1
    end
    file_out.close

  end


  def build_cluster_info clusters, genome_index

    puts "Building Cluster Info (Consensus).."
    # Clusters-Info.tsv
    file_out = File.open(@cdHitFile+".info.tsv","w")
    file_out.write("cluster_id\toccurences\tconsensus_length\tconsensus_name\tlengths_distribution\tnames_distribution\n")
    clusters.each do |c|
      lengths_count = c[:lengths].group_by{|x| x}.map{|k,v| [k, v.count] }
      gene_names_count = c[:gene_names].group_by{|y| y}.map{|k,v| [k, v.count] }
      lengths_sorted = lengths_count.sort_by { |a| a[1].to_i }.reverse
      gene_names_sorted = gene_names_count.sort_by { |a| a[1].to_i }.reverse
      file_out.write("cluster_#{c[:id]}\t#{c[:occurences]}/#{genome_index.size}\t#{lengths_sorted[0][0]}\t#{gene_names_sorted[0][0]}\t#{lengths_count.sort_by { |a| a[0].to_i }}\t#{gene_names_sorted.sort_by { |a| a[1].to_i }}\n")
    end
    file_out.close

  end

end




## Main ##

if __FILE__==$0

  # Define the usage of this script
  usage = "\n$ cd-hit-parser.rb <cdHitFile> <samplesFile>"
  usage += "\n\n  ! Results header should look like this to map the samples file"
  usage += "\n   >FM209186_PLES_28951"
  usage += "\n   >SeqID_LocusTag\n\n"

  if ARGV.length < 2
    abort usage
  end

  cdHitFile = ARGV[0]
  fastaFile = ARGV[1]
  cdhit_parser = CdhitParser.new(cdHitFile, fastaFile)
  genome_index, genome_genes = cdhit_parser.build_genome_index
  clusters = cdhit_parser.parse_cdhit_file
  cdhit_parser.build_cluster_matrix(clusters, genome_index, genome_genes)
  cdhit_parser.build_cluster_info(clusters, genome_index)

end




