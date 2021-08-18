#!/usr/bin/perl
# vim: syntax=perl tabstop=4 expandtab

#-----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 1, 2015
#-----------------------------------
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my ( $matrix, $row_names, $sample_names ) = get_matrix( $$options{ 'star_log' } );
print_matrix( $matrix, $row_names, $sample_names );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'star_log|f=s@', 'help|h' );
	unless( $$options{ 'star_log' } ) {
		my $usage = "$0 <--star_log|-f>";
		print STDERR $usage, "\n";
		exit 1;
	}
	return $options;
}

sub get_matrix {
	my( $log_files ) = @_;
	my $matrix = {};
	my $row_names = [];
	my $sample_names = [];
	foreach my $log_file( @$log_files ) {
		my( $sample_name ) = ( basename( $log_file ) =~ /(.+)\.Log\.final\.out/ ); 
		push @$sample_names, $sample_name;
		open( FH, "<$log_file" ) or die "Error in opening the file, $log_file, $!\n";
		my $flag = 0;
		while( my $line = <FH> ) {
			chomp $line;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			if( $line =~ /^Number/ && ! $flag ) {
				$flag = 1;
			}
			if( $flag ) {
				my( $key, $value ) = split(/\|/, $line );
				if( $value ) {
					$key =~ s/^\s+//;
					$key =~ s/\s+$//;
					$key =~ s/,//g;
					$key =~ s/\s+/_/g;
					$value =~ s/^\s+//;
					$value =~ s/\s+$//;
					push @$row_names, $key unless( exists $$matrix{ $key } );
					$$matrix{ $key }{ $sample_name } = $value;
				}
			}
		}
		close FH or die "Error in closing the file, $log_file, $!\n";
	}
	return ( $matrix, $row_names, $sample_names );
}

sub print_matrix {
	my( $matrix, $row_names, $sample_names ) = @_;
	# my $STAR_version = `STAR --version`;
	my $header = join( ",", ( "STAR", @$sample_names ) );
	print $header, "\n";
	foreach my $feature( @$row_names ) {
		my @values = ();
		FOR_LOOP:
		foreach my $sample( @$sample_names ) {
			unless( exists $$matrix{ $feature }{ $sample } ) {
				last FOR_LOOP;
			}
			push @values, $$matrix{ $feature }{ $sample } or undef;
		}
		unless( @values ) {
			print $feature, "\n";
		} else {
			print join( ",", ( $feature, @values ) ), "\n";
		}
	}
}