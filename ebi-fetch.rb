#!/usr/bin/env ruby
# ======================================================================
# WSDbfetch SOAP using soap4r 
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/dbfetch
# http://www.ebi.ac.uk/Tools/webservices/tutorials/ruby
# ======================================================================
# Load library
require 'bio'
require 'uri'
require 'open-uri'


usage = 
"
 ebi-fetch.rb 
    options =	db-list         		# list database available
		fetch-seq <db> <id> <out>	# fetch fasta sequence from db
		fetch-embl <db> <id> <out>	# fetch annotation from db
		fetch-struct <db> <id> <out>	# fetch structure from the pdb
"

opt = ARGV[0]
@server = Bio::Fetch.new('http://www.ebi.ac.uk/cgi-bin/dbfetch')

def listdb
  puts @server.databases
end

def fetch filter
  db = ARGV[1] or die usage
  id = ARGV[2] or die usage

  if ARGV[3]
    out = ARGV[3]
  else
    out = "#{id}"
  end

  if filter == 1
    f = File.new("#{out}.fasta", "w")
    f.write(@server.fetch(db,id,'raw','fasta'))
    f.close
  elsif filter == 2
    f = File.new("#{out}.embl", "w")
    f.write(@server.fetch(db,id,'raw','embl'))
    f.close
  elsif filter == 3
    f = File.new("#{out}.pdb", "w")
    f.write(@server.fetch(db,id,'raw','pdb'))
    f.close
  end
end



## Main ##

# switch case with the options
case opt

when "db-list"
  listdb
when "fetch-seq"
  fetch 1
when "fetch-embl"
  fetch 2
when "fetch-struct"
  fetch 3
else
  puts usage
end

