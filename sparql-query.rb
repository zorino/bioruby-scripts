#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	Maxime Déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	13-03-15
# version: 	0.01
# licence:  	

# Need a configuration file with the query

require 'rdf'
require 'sparql/client'

usage = "

sparql-query.rb <config file>

In config file : 	
	-endpoint 	set the endpoint URL after this tag
	-query		set the query afther this tag

Optional arguments :

	you can put any argument in the command line that you want
	to be replaced in your query

"

abort usage if ARGV.length == 0


def readArgument
  argList = {}
  for i in 1..(ARGV.length-1)    
    if ARGV[i].include? "="
      key = ARGV[i].split("=")
      argList[key[0].to_s] = key[1].to_s
    end
  end 
  return argList
end



def readConfig
  file = ARGV[0]
  argList = {}

  argList = readArgument

  fileInfo = File.open("sparql-result.info.txt", "a")

  f = File.open (file)

  inEndpoint = 0
  inQuery = 0
  @queryString = ""
  @endpointString = ""

  fileInfo.write("\n# # # Your Query # # #\n")

  f.each_line do |l|
    l.chomp!
    if l.empty?
      next
    elsif l.downcase.include? "-query"
      inQuery = 1
      inEndpoint = 0
      l.slice! "-query"
    elsif l.downcase.include? "-endpoint"
      inEndpoint = 1
      inQuery = 0
      l.slice! "-endpoint"
    end
    
    if inQuery==1 && l.strip != ""
      # if l.include? "SELECT"        
      # end
      argList.each do |k,v|
        # puts "key : #{k} \t value: #{v}"
        l.gsub!("#{k}","#{v}")
        fileInfo.write(l)
      end
      @queryString = @queryString + "\n" + l.strip
    end
    if inEndpoint==1
      @endpointString = @endpointString + l.strip
    end

  end

  f.close
  fileInfo.close

end


# Build solution by identifying each ARGMUMENT in the SELECT query
def buildSolution(symbol)

  solArg = []
  solOutput = ""
  queryArray = @queryString.split("\n")

  queryArray.each do |l|
    if l.downcase.match(/^select/)
      solArg = l.scan(/\w+/)
    end
  end

  for i in 1..solArg.length-1
    if solOutput != ""
      solOutput = solOutput + "+ \"\\t\" + " + " #{symbol}[:#{solArg[i]}]"
    else
      solOutput = "#{symbol}[:#{solArg[i]}]"
    end
  end

  return solOutput

end                             # end of buildSolution


def sendSparqlQuery

  fileInfo = File.open("sparql-result.info.txt", "a")

  fileInfo.write("\n# # # RESULTS # # #\n")
  
  solOutput = buildSolution("solution")

  endpoint = SPARQL::Client.new(@endpointString)

  endpoint.query(@queryString+'&format=xml').each do |solution|

    # if CONSTRUCT only
    fileInfo.write(solution.to_s)

    #  if SELECT only
    out = ""
    solution.each_binding do |n,v|
      out << "#{v}\t"
    end
    puts out
    fileInfo.write(out)

  end

  fileInfo.close

end


# main #

readConfig
sendSparqlQuery

