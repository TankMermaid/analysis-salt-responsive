#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'arginine'
require 'set'

par = Arginine::parse do
  desc ""
  arg :pval
  arg :fc
end

def file2hash(fn)
  File.readlines(fn).map do |line|
    otu, val = line.chomp.split("\t")
    [otu, val.to_f]
  end.drop(1).to_h
end

pval_hash = file2hash(par[:pval])
fc_hash = file2hash(par[:fc])

raise unless pval_hash.keys.sort == fc_hash.keys.sort
otus = pval_hash.keys.sort
pvals = pval_hash.values_at(*otus)
fcs = fc_hash.values_at(*otus)

max_pval = 0.05
min_opacity = 0.0
target_size = 50.0
max_fc = fcs.map(&:abs).max

puts "DATASET_SYMBOL"
puts "SEPARATOR SPACE"
puts "DATASET_LABEL foo"
puts "COLOR #ff0000"
puts "MAXIMUM_SIZE #{target_size}"
puts "DATA"

otus.zip(pvals, fcs).each do |otu, pval, fc|
  min_log = 3
  opacity = 1.0 - (Math.log10(pval) + min_log) / min_log

  if pval < max_pval
    fill = 1
  else
    fill = 0
    # swo> make lines all dark
    opacity = 1.0
  end

  if fc < 0
    rgb = "255,0,0"
  else
    rgb = "0,0,255"
  end

  color = "rgba(#{rgb},#{opacity})"
  size = target_size * (fc.abs / max_fc)
  puts ["'" + otu + "'", 2, size, color, fill, 1.0].join(" ")
end
