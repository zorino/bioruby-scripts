#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	13-06-10
# version: 	0.01
# licence:  	

# Dependencies : libxml-ruby

require 'bio'

# Convert newick to phyloxml with forester.jar (https://code.google.com/p/forester/)
# 
# java -cp path\to\forester.jar org.forester.application.phyloxml_converter -f=<field option> [options] <infile> <outfile>
# see http://goo.gl/UrXnmC


class BioTreeManip

  attr_accessor :tree

  def initialize tree_file
    @tree = Bio::PhyloXML::Parser.open(tree_file).next_tree
  end


  def node_iterator
    @tree.each_node do |node|
      p node.name.to_s + " " + node.confidences.to_s + " " + node.other.to_s + "\n"
    end
  end

  # annotation_hash[sample] = {'desc': "description string", 'uri': "http//genome.url"}
  def add_annotation_on_nodes annotation_hash

    @tree.each_node do |node|

      if annotation_hash.has_key? node.name

        annotation = Bio::PhyloXML::Other.new()
        annotation.element_name = "annotation"

        annotation_hash[node.name].each do |k,v|
          description = Bio::PhyloXML::Other.new()
          description.element_name = "#{k}"
          description.value = "#{v}"
          annotation.children  << description
        end

        node.properties << annotation

      end

    end

  end


  # property_hash[sample] = {'attributes': {'ref'=>'style:node_color'}, 'value': "#ff00cc"}
  # <property ref="style:node_color" datatype="xsd:token" applies_to="node">#ff00cc</property>
  def add_property_on_nodes property_hash

    @tree.each_node do |node|

      if property_hash.has_key? node.name

        property_hash[node.name].each do |a| 
          property = Bio::PhyloXML::Other.new()
          property.element_name = "property"
          property.value = a['value']
          property.attributes = a['attributes']
          node.properties << property
        end

      end

    end

  end


  # chart_hash[sample] = {'chart_type-1' => "value(color)", 'chart_type-2' => "value(color)"}
  def add_chart_on_nodes chart_hash

    @tree.each_node do |node|

      if chart_hash.has_key? node.name

        chart = Bio::PhyloXML::Other.new()
        chart.element_name = "chart"

        chart_hash[node.name].each do |k,v|
          chart_type = Bio::PhyloXML::Other.new()
          chart_type.element_name = "#{k}"
          chart_type.value = "#{v}"
          chart.children  << chart_type
        end

        node.properties << chart

      end

    end

  end

  # change name with possible attributes
  def change_node_name *attribute

    @tree.each_node do |node|
      next if node.name.nil?
      sample_name = node.name.to_s.dup
      node.name = nil
      new_name = Bio::PhyloXML::Other.new()
      new_name.element_name = "name"
      if attribute[0]
        if attribute[0].has_key? sample_name
          # new_name.attributes = {'bgStyle' => "bold"}
          attribute_hash = {}
          attribute[0][sample_name].each do |k,v|
            attribute_hash[k] = v
          end
          new_name.attributes = attribute_hash
        end
      end
      new_name.value = sample_name
      node.properties << new_name
    end

  end

  # Write phyloxml tree to file
  def write_tree outfile
    writer = Bio::PhyloXML::Writer.new(outfile)
    writer.write(@tree)
  end

  # Write the phyloxml file into a js var object
  def print_json_obj
    
  end


end

###---------- END of class


# EXAMPLE

# def create_annotation_hash file

#   hash = Hash.new

#   File.open(file,"r") do |fileR|
#     while l = fileR.gets
#       lA = l.chomp!.split("\t")
#       sample = lA[0]

#       # Sample
#       description =  "Sample:    #{sample}   #{lA[26]}\\n"

#       # Patient
#       sex = "Male"
#       if lA[25].to_i == 2
#         sex = "Female"
#       end
#       description += "Patient:    #{sex} #{lA[24]} yold\\n"

#       # Sample MLST
#       description += "MLST:     ST-#{lA[27]}\\n"

#       # Toxin
#       description += "Toxins:     "
#       description += "ToxinA," if lA[2].include? "pos"
#       description += "Toxin1," if lA[6].include? "pos"
#       description += "Tpi," if lA[3].include? "pos"
#       description += "CdtA," if lA[4].include? "pos"
#       description += "CdtB," if lA[5].include? "pos"
#       description += "TcdA," if lA[8].include? "pos"
#       description += "TcdB," if lA[9].include? "pos"
#       description += "TcdC," if lA[7].include? "pos"
#       if description[-1] == ","
#         description[-1] = ""
#       end
#       description += "\\n"

#       # Resistance
#       description += "Resistance:  "

#       description += "Azi," if lA[10].to_i > 8
#       description += "Cax," if lA[11].to_i > 32
#       description += "Cla," if lA[12].to_i > 4
#       description += "Cd," if lA[14].to_i > 4
#       description += "MNZ," if lA[17].to_i > 2
#       description += "Mox," if lA[18].to_i > 4
#       description += "Va," if lA[20].to_i > 2
#       if description[-1] == ","
#         description[-1] = ""
#       end
#       description += "\\n"

#       # Assembly
#       description += "Assembly:  #{lA[21]} nt in #{lA[22]} contigs"
#       description += "\\n"


#       hash[lA[0]] = {'desc' => description}

#     end
#   end

#   hash

# end



# # Main #
# tree = BioTreeManip.new(ARGV[0])

# ## get annotation from metadata
# annotation_hash = create_annotation_hash ARGV[1]
# tree.add_annotation_on_nodes annotation_hash

# # create barChart for virulence resistance based on metadata
# chart_hash = {}
# annotation_hash.each_key do |k|
#   desc_array = annotation_hash[k]['desc'].split("\\n")
#   if desc_array[3].to_s.include? ","
#     virulence = (desc_array[3].to_s.split(",").length.to_f*100/8).round(2)
#   else
#     virulence = 0
#   end
#   if desc_array[4].to_s.include? ","
#     resistance = (desc_array[4].to_s.split(",").length.to_f*100/7).round(2)
#   else
#     resistance = 0
#   end

#   mlst = "st_" + desc_array[2].to_s.gsub("MLST:     ST-","").chomp
#   chart_hash[k] = {'mlst' => mlst.to_s}

#   chart_hash[k] = {'virulence' => virulence.to_s,
#                    'resistance' => resistance.to_s,
#                    'mlst' => mlst}
# end
# tree.add_chart_on_nodes chart_hash 

# # # Add MLST bin chart
# # stchart_hash = {}
# # annotation_hash.each_key do |k|
# #   desc_array = annotation_hash[k]['desc'].split("\\n")
# #   mlst = "st_" + desc_array[2].to_s.gsub("MLST:     ST-","").chomp
# #   chart_hash[k] = {'mlst' => mlst.to_s}
# # end
# # tree.add_chart_on_nodes stchart_hash 

# # Change Node Names
# node_names = {}
# tree.tree.each_node do |n|
#   if ! n.name.nil?
#     node_names[n.name] = {'bgStyle' => n.name.split("-")[0]}
#   end
# end

# tree.change_node_name node_names
# tree.write_tree "ztree.xml"
