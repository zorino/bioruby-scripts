#!/usr/bin/env ruby

require 'bio'

# SOME EXAMPLES :
# ncbi = Bio::NCBI::REST.new
# ncbi.efetch("185041", {"db"=>"nucleotide", "rettype"=>"gb", "retmode" => "xml"})
# ncbi.efetch("J00231", {"db"=>"nuccore", "rettype"=>"gb", "retmode"=>"xml"})
# genbank = ncbi.efetch("AAA52805", {"db"=>"protein", "rettype"=>"gb"})

# Set argumets

usage = 
"
Usage : ncbi-fetch.rb <database> <id> <output-type> [output name]
    options :	info		gives databases information
		database  	nt or prot
		output type = 	fasta -> seq only
				gb -> genbank
"

def run_ncbi_fetch db, id, type, out

  Bio::NCBI.default_email = 'default@default.com'
  ncbi = Bio::NCBI::REST.new

  ids = id.split(",")
  db = "protein" if db.include? "prot"
  db = "nucleotide" if db.include? "nt"

  tryout = 0

  begin
    genbank = ncbi.efetch(ids, {"db"=>db, "rettype"=>type})
    if genbank.include? "Cannot process ID list"
      abort "#{id} doesn't seem to exist in <#{db}> database"
    elsif genbank.include? "<!DOCTYPE html"
      puts "server request failed"
      fail
    else
      f = File.new(out, "w")
      f.print(genbank)
      f.close
    end
  rescue
    if tryout < 5
      tryout += 1
      sleep (tryout*5)
      retry
    end
  end

  return genbank

end


if ARGV[0] == "info"
  puts "== Databases Available =="
  puts ncbi.einfo()
  puts "== =="
  abort()
end

db = ARGV[0] or abort usage
id = ARGV[1] or abort usage
type = ARGV[2] or abort usage
out = ""
if ! ARGV[3].nil?
  out = ARGV[3]
elsif type == "gb"
  out = "#{id}.gbk"
else
  out = "#{id}.#{type}"
end

puts "# Fetching #{id}"
genbank = run_ncbi_fetch db, id, type, out


if type == "gb" and (genbank =~ /^WGS/)
  Dir.mkdir(id) if ! Dir.exists? id
  genbank.split("\n").each do |l|
    if l[0..2] == "WGS"
      puts "# This is a WGS ! ..fetching sub genbank !"
      ids = l.split(/\s+/)[1].split("-")
      if ids.length == 1
        run_ncbi_fetch db, ids[0], type, "#{id}/#{ids[0]}.gbk"
      else
        for i in ids[0][4..-1]..ids[1][4..-1]
          tmp_id = ids[0][0..3]+i.to_s
          puts "  ..fetching #{tmp_id}"
          run_ncbi_fetch db, tmp_id, type, "#{id}/#{tmp_id}.gbk"
        end
      end
    end
  end
end
