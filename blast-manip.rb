#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-01-10
# version: 	0.01
# licence:  	

require 'bio'


usage = 
"
blast-manip.rb id=0.9 [ ]  <blastOutput> 
	Options		covLength	# lenght of coverage
			covQ		# % coverage on Query
			covS		# % coverage on Subject
			besthit		# return the 1st (best) hit

"

abort usage if ARGV.length == 0

report = Bio::FlatFile.auto(ARGV[ARGV.length-1]) or abort usage

idOpt = ARGV[0].split("=") or abort usage
if idOpt[1].to_f <= 1
  id = idOpt[1].to_f*100
else
  id = idOpt[1].to_f
end
puts "ID opt = #{id}"


# Iterate through Query->Hit->HSPs
queryResults = Hash.new{|h,k| h[k] = []}
hitsArray = []
hitPos = {}
uniqueRegion = 1                # HSP is not included in a previous one

report.each_entry do |entry|
  queryName = entry.query_def
  entry.each_hit do |hit|
    hitPos.clear
    align_length = 0
    hitName = hit.definition
    hit.each do |hsp|
      if hsp.percent_identity.to_f >= id
        if ! hitPos.empty?
          hitPos.each_key do |k|
            if (hsp.query_from < hitPos[k] and hsp.query_to > k)
              uniqueRegion = 0
            end
          end
          hitPos[hsp.query_from] = hsp.query_to if uniqueRegion == 1
        else
          hitPos[hsp.query_from] = hsp.query_to
        end
        align_length += hsp.align_len if uniqueRegion == 1
        uniqueRegion = 1
      end
    end                         # end of hsp
    hitsArray.push("#{hitName}\t#{align_length}") if align_length > 0
    align_length = 0
  end                           # end of hit
  queryResults[queryName].push(*hitsArray)
  hitsArray.clear
end                             # end of query

# Print the result
queryResults.each_key do |k|
  queryResults[k].each do |v|
    puts "#{k}\t#{v}"
  end
end

report.close


