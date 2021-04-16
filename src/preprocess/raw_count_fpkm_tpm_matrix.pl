#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my $meta_info = exists $$options{ 'metasheet' } ? get_metasheet_info( $$options{ 'metasheet' } ) : undef;
my ( $matrix, $files ) = get_matrix( $options );
print_matrix( $matrix, $options, $files, $meta_info );
exit $?;


sub parse_options {
	my $options = {};
	GetOptions( $options, 'file|f=s@', 'column|l:i', 'htseq|t', 'cufflinks|c', 'salmon|s','remove_ERCC_ids|e', 'header|d', 'metasheet|m:s', 'help|h' );
	unless( $$options{ 'file' } ) {
		my $usage = "$0 <file|-f> [--file|-f] [--column|-l <1 or 2 or 3; default=3>] [--htseq|-t] 
				[--header|-d <having this option removes header line>] [--cufflinks|-c] [--remove_ERCC_ids|-e]
                [--metasheet|-m <having this option generates a matrix with header from metasheet annotation columns";
		print STDERR $usage, "\n";
		exit 1;
	}
	unless( $$options{ 'column' } ) {
		$$options{ 'column' } = 3;
	}
	return $options;
}

sub get_metasheet_info {
    my( $file ) = @_;
    my $info = {};
    open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
    my @columns = ();
    my $header = <FH>;
    chomp $header;
    my @head = split(",", $header);
    foreach my $index( 1 .. scalar @head -1 ) {
        unless( substr($head[$index], 0, 5) eq 'comp_' ) {
            push @columns, $index;
        }
    }
    my $count = 1; 
    while( my $line = <FH> ) {
        chomp $line;
        my @array = split(",", $line );
        foreach my $col( @columns ) {
            $array[$col] =~ s/_/./g;
        }
        $$info{ $array[0] } = join("_", @array[@columns]) . '_' . $count;
        $count++;
    }
    close FH or die "Error closign the file, $file, $!\n";
    return $info;
}

sub get_matrix {
	my( $options ) = @_;
	my $basenames = [];
	my $matrix = {};
	my @files = @{ $$options{ 'file' } };
	foreach my $file( @files ) {
		open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
		my $file_base  = basename( $file );
		if( $$options{ 'cufflinks' } ) {
			$file_base = basename( dirname( $file ) );
		} elsif( $$options{ 'htseq' } ) {
			$file_base =~ s/\.htseq\.counts//;
		} elsif( $$options{ 'salmon' } ){
			$file_base =~ s/\.quant\.sf//;
		} else {
			$file_base =~ s/\.counts\.tab//; #STAR output
			#------ RSEM out files contain these ----------#
			$file_base =~ s/\.genes\.results//;
			$file_base =~ s/\.isoforms\.results//;
			#----------------------------------------------#
		}
		push @$basenames, $file_base;
		my $header = <FH> if $$options{ 'header' };
		while( my $line = <FH> ) {
			chomp $line;
			if( $$options{ 'cufflinks' } ) {
				my @array_of_vals = split( "\t", $line );
				#PARSE cufflinks output- take gene_id and FPKM
				my( $gene_id, $fpkm ) = ( $array_of_vals[ 3 ], $array_of_vals[ 9 ] );
				if( not substr( $gene_id, 0, 2 ) eq '__' ) {
					if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
						# do nothing
					} else {
						$$matrix{ $gene_id }{ $file_base } = $fpkm;
					}
				}
			} elsif($$options{ 'salmon' }) {
				my @array_of_vals = split( "\t", $line );
				#PARSE cufflinks output- take gene_id and TPM
				my( $gene_id, $tpm ) = ( $array_of_vals[ 0 ], $array_of_vals[ 3 ] );
				if( not substr( $gene_id, 0, 2 ) eq '__' ) {
                                        if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
                                                # do nothing
                                        } else {
						$$matrix{ $gene_id }{ $file_base } = $tpm;
					}
				}
			} elsif( $$options{ 'htseq' } ) {
				my( $gene_id, $count ) = split( "\t", $line );
				if( not substr( $gene_id, 0, 2 ) eq '__' ) {
                                        if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
                                                # do nothing
                                        } else {
						$$matrix{ $gene_id }{ $file_base } = $count;
					}
				}
			} else {
				my @array_of_vals = split( "\t", $line );
				my( $gene_id, $count ) = @array_of_vals[ 0, $$options{ 'column' } ];
                                if( not substr( $gene_id, 0, 2 ) eq 'N_' ) {
                                        if( $$options{ 'remove_ERCC_ids' } && $gene_id =~ /ERCC\-00\d\d\d/ ) {
                                                # do nothing
                                        } else {
                                                $$matrix{ $gene_id }{ $file_base } = $count;
                                        }
                                }

			}
		}
		close FH or die "Error in closing the file, $file, $!\n";
	}
	return $matrix, $basenames;
}

sub print_matrix {
	my( $matrix, $options, $files, $meta_info ) = @_;
	my @header = $meta_info ? get_seurat_header($files, $meta_info) : @$files;
    print STDOUT join(",", ("Gene_ID",@header)), "\n";
	foreach my $gene_id( keys %$matrix ) {
		my @counts = ();
		foreach my $file_base( @$files ) { 
			if( exists $$matrix{ $gene_id }{ $file_base } ) {
				push @counts, $$matrix{ $gene_id }{ $file_base };
			}
			else {
				push @counts, 0;
			}
		}
		print STDOUT join( ",", ( $gene_id, @counts ) ), "\n";
	}
}

sub get_seurat_header {
    my( $files, $meta_info ) = @_;
    my @header = ();
    foreach my $file( @$files ) {
	#BUG FIX: filenames are NOT sample names-
	#need to parse out sample names
	my $tmp = (split(/\ /, $file))[0];
        push @header, $tmp;
    }
    return @header;
}
