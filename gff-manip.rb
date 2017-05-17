#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-07-26
# version: 	0.01
# licence:  	

require 'bio'


class GFF2GBK

  # Initialize the GFF2GBK object with a list of gff files.
  # Load everything in gffFilesHash (concurrently).
  def initialize gffFilesToProcess, directoryIn
    
    @directory = directoryIn
    maxThreads = 64
    @gffFilesHash = Hash.new
    myThreads = Array.new
    
    # while processedKeys.size < gffFilesToProcess.size
    while gffFilesToProcess.size > 0
      
      if gffFilesToProcess.size < maxThreads
        maxThreads = gffFilesToProcess.size
      end
 
      maxThreads.times do |i|

        myThreads[i] = Thread.new do
          
          curKey = gffFilesToProcess.shift
          Thread.current[:name] = curKey.to_s
          
          @gffFilesHash[curKey] = Bio::GFF::GFF3.new(File.open("#{directoryIn}/#{curKey}.gff"))
          
        end

      end
          
      # myThreads.each { |t| t.join ; puts "Thread #{t[:name]} finished !" }             
      myThreads.each { |t| t.join }             
        
    end  

  end                           # Initialize
  

  
  # Create a Hash of all Genbank files parsed from gff files
  def createGbkFiles
    
    @gbkFilesHash = Hash.new
    
    @gffFilesHash.each do |key, gff|
     
      gbkString = createOneGbk(key,gff)
      @gbkFilesHash[key] = gbkString
      
      # puts @gbkFilesHash[key]
      # puts "Create #{key} file !"
      
    end
    
  end                           # Method : createGbkFiles
  


  # Dumb all the Gbk files when created if not will create it.
  # @gbkFilesHash
  def dumpGbkToDisk directoryOut
    
    # puts "Dumping all gbk files pls wait ..."

    Dir.mkdir(directoryOut) unless Dir.exist? directoryOut
    Dir.chdir(directoryOut)

    @gbkFilesHash.each do |key,gbk|

      File.open("#{key}.gbk","w") do |fopen|
        fopen << gbk
      end

    end

  end                           # Method : dumpAllGbk


  # # #
  # Private Methods
  # # #

  # private
  def grepFastaSeq file
  
    fasta = ""

    File.open(file) do |f|
      state = 0
      while l = f.gets
        
        fasta += l.to_s unless state == 0
        
        if l[0] == ">"
          state = 1
        end     
        
      end
      
    end
    fasta
    
  end                           # Private Method : grepFastaSeq

  
  # private
  def createOneGbk key, gff
    
    featuresArray = []
    seqname = ""
    fastaSeq = grepFastaSeq "#{@directory}/#{key}.gff"
    
    # puts gff.to_s

    gff.records.each do |rec|

      prod = ""
      loc = ""

      if seqname == "" and rec.seqname   
        seqname = rec.seqname 
      end
      
      if rec.feature == "CDS"
        
        loc = "#{rec.start}..#{rec.end}"

        if rec.strand == "-"
          loc = "complement(#{loc})"
        end

        # puts "Location : #{loc}"
        
        qualifiersArray = []
        
        rec.attributes.each do |att|
          
          # qualif = k.split("=")
          
          # i = 1
          # while i < qualif.size
          #   qualifString += qualif[i].to_s
          # end

          # if rec.attributes[k]
          #   qualifString = qualifString + " " + rec.attributes[k]        
          # end       

          # qualifiersArray.push(Bio::Feature::Qualifier.new(qualif[0].to_s,qualifString))

          qualifiersArray.push(Bio::Feature::Qualifier.new(att[0],att[1]))         
          
        end

        newFt = Bio::Feature.new("CDS", loc ,qualifiersArray)

        featuresArray.push(newFt)
      end
    end

    newBioSeq = Bio::Sequence.new(fastaSeq.to_s)
    newBioSeq.entry_id = seqname
    newBioSeq.features = featuresArray

    # puts newBioSeq.output(:fasta)

    newGbk = Bio::GenBank.new(newBioSeq.to_s)

    gbkString = newBioSeq.output(:GenBank).to_s.force_encoding("UTF-8").gsub /^$\n/, ''

    gbkString # return Genbank file as string
  
  end                           # Private Method createOneGbk
  

  # declare private methods
  private :grepFastaSeq, :createOneGbk


end                             # CLASS


# Main #

usage = "

  $ gff-manip.rb <gff file>

"

directory = ARGV[0]

gffFilesToProcess = Array.new
Dir.foreach(directory) do |f|
  next if f == '.' or f == ".."
  gffFileName = f.gsub!(".gff","")
  gffFilesToProcess.push gffFileName
end


puts "Initializing GFF files ..."
gff2gbk = GFF2GBK.new(gffFilesToProcess,directory)
puts "Creating Gbk Files"
gff2gbk.createGbkFiles
puts "Dumping All Gbk files to disk..."
gff2gbk.dumpGbkToDisk "#{directory}-Gbk"


# # # # # 
# redis #
# # # # #

# # #
#  Fct : while parsing dump on a redis.io in memory
#
def dump2redisio


end

# # #
# Fct : load from redis r/w
#
def write2redisio


end
