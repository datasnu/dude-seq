#!/usr/bin/perl -w

use POSIX;
use strict;
use warnings;

my $flow_order = "TACG";

while (my $id_line = <STDIN>)
{
	my $int_line = <STDIN>;
	
	chomp($id_line);
	chomp($int_line);
	
	if ($int_line =~ /^(\d+)\s+(\d+)\s+(.+)$/)
	{
		my $end_pos = $2;
		my @flow_arr = split(/\s+/, $3);
		my $called_seq = BaseCall($end_pos, @flow_arr);
		$called_seq =~ s/^TCAG//;
		
		print "$id_line\n";
		print "$called_seq\n";
	}
}


sub BaseCall
{
	my ($end_pos, @flowgram) = @_;
	my $called_seq = "";
	
	for (my $i = 0; $i < $end_pos; $i++)
	{
		my $signal = floor($flowgram[$i] + 0.5);
		my $base = substr($flow_order, $i % 4, 1);
		
		for (my $j = 0; $j < $signal; $j++)
		{
			$called_seq .= $base;
		}
	}
	
	return $called_seq;
}
