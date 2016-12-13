#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'arginine'

par = Arginine::parse do
  desc ""
  arg :taxa
  opt :template, default: "dataset_color_strip_template.txt"
end

template = File.open(par[:template]).read
template.gsub!("#SHOW_INTERNAL 0", "SHOW_INTERNAL 1")

taxa = File.readlines(par[:taxa]).map do |line|
  "'" + line.chomp + "'"
end

puts template

taxa.each do |taxon|
  puts [taxon, "#000000", taxon].join(" ")
end
