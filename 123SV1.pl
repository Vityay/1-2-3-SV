#!/usr/bin/perl -w
use strict;
use Config::Std;

# 1-2-3-SV step 1
# v 0.04
# Last modified 20120506

my %cfg = ();
die "Cannot read configuration file 123SV.conf" unless -r '123SV.conf';
read_config( '123SV.conf' => %cfg );
check_parameters();

my %lib_ids = ();
while ( my $lib = shift ) {
    die "Wrong library number given: $lib" unless $lib =~ m/^(\d+)$/;
    warn "Warning: library $1 was mentioned twice, a typo?\n" if exists( $lib_ids{$1} ); 
    $lib_ids{$1} = 1; 
}

# Determine which libraries to preocess
unless ( keys %lib_ids ) {
    my $answer = q{};
    while ( $answer !~ m/^[yn]/) {
        warn "No library id is given, process all? [y/n]?\n";
        $answer = <STDIN>;
        chomp $answer;
        if ( $answer eq 'y' ) {
            foreach my $lib_id ( keys %cfg ) {
                next unless $lib_id =~ m/^\d+$/;
                if (-e 'sv_size_distribution.'.$cfg{$lib_id}{'name'} ) {
                     warn "lib $lib_id seems to be processed, skipping it\n";
                     next;
                }
                $lib_ids{$lib_id} = 1 if $lib_id =~ m/^(\d+)$/; 
            }
        }
    }
}

my %seen  = (); # Clonality filtering 
my %fwd   = (); 
my %rev   = (); 
my %distr = (); 
my %stats = ();
my $lines_read = 0;
my $lib_id = q{};

foreach $lib_id ( sort {$a<=>$b} keys %lib_ids ) { 

    if (-e 'sv_size_distribution.'.$cfg{$lib_id}{'name'} ) {
         warn "lib $lib_id seems to be processed, skipping it\n";
         next;
    }

    die 'Incorrect number of chunks is given for library ',$lib_id, "\n" if $cfg{$lib_id}{'chunks'} !~ m/^[1-9]\d*$/;

    undef %distr;
    undef %stats;
    undef %seen;
    undef %fwd;
    undef %rev;
    
    read_paired( $lib_id );
    save_distr(  $lib_id );
    save_stats(  $lib_id );
}


sub check_parameters {
    die "Cannot read General::max_distance from configuration file" unless exists( $cfg{'Distribution'}{'max_distance'} );
    die "Define parameter remove_clonal in section Distribution" unless exists( $cfg{'Distribution'}{'remove_clonal'} );
}

sub read_paired {

    my $lib_id = shift;
    my $lines_read = 0;

    my $lib = $cfg{$lib_id}{'name'};
    my $libprefix = $cfg{$lib_id}{'prefix'};
    my $lib_type  = uc( $cfg{$lib_id}{'type'} ) or die "Library type for library $lib_id not specified";
    warn 'Determining insert sizes for lib # ', $lib_id, ' : type ', $lib_type, ' in ', $cfg{$lib_id}{'chunks'}, " chunks\n";
    open FDIST, '>', 'sv_size_distribution.'.$lib;
    open FSTAT, '>', 'sv_pairs_stats.'.$lib;
    open FNUCR, '>', 'sv_pairs_remote.'.$lib;
    open FNUCI, '>', 'sv_pairs_inversions.'.$lib;
    open FNUCE, '>', 'sv_pairs_evertions.'.$lib;
    foreach my $chunk ( 0 .. ( $cfg{$lib_id}{'chunks'} - 1 ) ) {
        undef %fwd;
        undef %rev;
        warn 'Started with chunk ', ($chunk + 1), ' of ', $cfg{$lib_id}{'chunks'}, ' for lib # ', $lib_id,' : ', $lib,' type ', $lib_type, "\n";
        open F1, $cfg{'General'}{'samtools_exe'}.' view '.$cfg{$lib_id}{'fileF'}.'|' or die 'Could not open file with forward tags for lib '.$lib;
        open F2, $cfg{'General'}{'samtools_exe'}.' view '.$cfg{$lib_id}{'fileR'}.'|' or die 'Could not open file with reverse tags for lib '.$lib;
        while ( !eof(F1) or !eof(F2) ) {
            # Forward read
            unless (eof(F1)) {
                my $line1 = <F1>;
                my ( $chr, $pos, $strand, $clone ) = get_coord( $lib_id, $line1 );
                report_stats( $chr, $pos, $lib_id, $chunk, $cfg{$lib_id}{'chunks'} ) if $chr and not(++$lines_read % 1_000_000);
                
                # Check whether we need this clone 
                my $in_chunk = 1;
                if ( $cfg{$lib_id}{'chunks'} > 1 ) {
                    my $sum = 0;
                    map { $sum += ord($_) } split //, $clone;
                    $in_chunk = 0 unless $sum % $cfg{$lib_id}{'chunks'} == $chunk;
                    
                }

                if ( $chr and $in_chunk ) { # Unambig mapped tag
                    die "error: multiple mapping position of forward tag for clone $clone" if exists( $fwd{ $clone } );

                    $fwd{$clone} = join("\t", $chr, $strand, $pos ); 
                    process_ditag( $clone, $lib_type ) if ( exists( $rev{$clone} ) );
                }
            }

            # Reverse read
            unless ( eof(F2) ) {
                my $line1 = <F2>;
                my ( $chr, $pos, $strand, $clone ) = get_coord( $lib_id, $line1 );
                report_stats( $chr, $pos, $lib_id, $chunk, $cfg{$lib_id}{'chunks'} ) if $chr and not( ++$lines_read % 1_000_000);
                
                # Check whether we need this clone 
                my $in_chunk = 1;
                if ( $cfg{$lib_id}{'chunks'} > 1 ) {
                    my $sum = 0;
                    map { $sum += ord($_) } split //, $clone;
                    $in_chunk = 0 unless $sum % $cfg{$lib_id}{'chunks'} == $chunk;
                    
                }

                if ( $chr and $in_chunk ) { # Unambig mapped tag
                    die "error: multiple mapping position of forward tag for clone $clone" if exists( $rev{ $clone } );

                    $rev{$clone} = join("\t", $chr, $strand, $pos );
                    process_ditag( $clone, $lib_type ) if exists($fwd{$clone});
                }
            }
        } # while !eof(F1) or !eof(F2)
        close F1;
        close F2;
        warn "Finished with chunk ",$chunk+1," of ".$cfg{$lib_id}{'chunks'}." for lib # $lib_id : $lib, type $lib_type\n";
    } # foreach chunk
}

sub get_coord {
    my ( $lib, $line ) = @_;
    my ( $clone, $map_flag, $c, $p, $qual, $match ) = split /\t/, $line;

    # Skip if prefix unless multi-lib BAM and prefix is specified
    if ( $cfg{$lib}{'prefix'} ) {
         my $prefix = $cfg{$lib}{'prefix'};
         return unless $clone =~ m/$prefix/;
    }

    if ( $map_flag & 4 ) {
        delete $fwd{$clone} if exists($fwd{$clone});
        delete $rev{$clone} if exists($rev{$clone});
    }

    return if exists($cfg{'Distribution'}{'max_ambiguity'} ) and $line =~ m/X0\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'max_ambiguity' };
    return if exists($cfg{'Distribution'}{'secondary_hits'}) and $line =~ m/X1\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'secondary_hits'};
    return if exists($cfg{'Distribution'}{'alignment_gaps'}) and $line =~ m/XO\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'alignment_gaps'};
    return if exists($cfg{'Distribution'}{'alignment_mismatches'}) and $line =~ m/XM\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'alignment_mismatches'};

    my $strand = '+';
    if ( $map_flag & 16 ) {
        die "Unexpected match pattern: $match" unless $match =~ m/^[\dDIMS]+$/;
        my $shift = 0;
        while ( $match =~ m/(\d+)([DIMS])/g ) {
            $shift += $1 if $2 eq 'M' or $2 eq 'D';
        }
        $p += $shift - 1;
        $strand = '-';
    }
    # Trim when needed
    foreach my $trim ( $cfg{$lib}{'trim'} ) {
        $clone =~ s/$trim// if $trim;
    }
    return ( $c, $p, $strand, $clone );
}

sub process_ditag {
    my $clo = shift;
    my $lib_tp = shift;
    my ( $c1, $s1, $p1 ) = split /\t/, $fwd{$clo};
    my ( $c2, $s2, $p2 ) = split /\t/, $rev{$clo};

    $stats{'total'}++;
    return if $cfg{'Distribution'}{'remove_clonal'} and ++$seen{ join( "\t", $c1, $s1, $p1, $c2, $s2, $p2 )} > 1; # we do not take it into account 
    if ( $c1 ne $c2 or abs($p1 - $p2 ) > $cfg{'Distribution'}{'max_distance'} ) { # remotely located di-tags
        $stats{'remote'}++;
        $stats{'nonclonal'}++;
        ( $c1, $c2, $p1, $p2 ) = ( $c2, $c1, $p2, $p1 ) if $c1 lt $c2 or ( $c1 eq $c2 and $p1 > $p2 );
        print FNUCR join( "\t", $c1, int($p1/$cfg{'General'}{'bin_size'}), $c2, int($p2/$cfg{'General'}{'bin_size'})), "\n";
    } 
    else { # local
        my $ori = ditag_ori($p1, $p2, $s1.$s2, $lib_tp);
        if ( $p2 < $p1 ) {
            $ori = reverse(lc($ori));
        }    
        $stats{$ori}++;
        $stats{'local'}++; 
        if (uc($ori) eq 'TH') { # correctly oriented
            $distr{abs($p1-$p2)+1}++;
            $stats{'consistent'}++;
            $stats{'balance'}++ if $ori eq 'TH';
            $stats{'balance'}-- if $ori eq 'th';
             
        }
        elsif ( uc($ori) eq 'HT' ) { # evertion
            ( $p1, $p2 ) = ( $p2, $p1 ) if $p1 > $p2; 
            print FNUCE join( "\t", $c1, int($p1/$cfg{'General'}{'bin_size'}), int($p2/$cfg{'General'}{'bin_size'})), "\n";
        }
        else { #inversion
            ( $p1, $p2 ) = ( $p2, $p1 ) if $p1 > $p2; 
            print FNUCI join( "\t", $c1, int($p1/$cfg{'General'}{'bin_size'}), int($p2/$cfg{'General'}{'bin_size'})), "\n";
        }     
        $stats{'nonclonal'}++;
    }
    delete $fwd{$clo};
    delete $rev{$clo};
}

sub ditag_ori {
    my ( $pos1, $pos2, $s12, $type ) = @_;
    my $class = q{};
    if ( $type eq 'PE' ) {
        $class = $s12 eq '+-' ? 'TH' : $s12 eq '++' ? 'TT' : $s12 eq '-+' ? 'HT' : $s12 eq '--' ? 'HH' : '??'; 
    }
    elsif ( $type eq 'MP_SOLID' ) {
        $class = $s12 eq '+-' ? 'HH' : $s12 eq '++' ? 'HT' : $s12 eq '-+' ? 'TT' : $s12 eq '--' ? 'TH' : '??'; 
    }
    elsif ( $type eq 'MP_SOLEXA' ) {
        $class = $s12 eq '+-' ? 'HT' : $s12 eq '++' ? 'HH' : $s12 eq '-+' ? 'TH' : $s12 eq '--' ? 'TT' : '??'; 
    }
    else {
        die "Wrong lib type: $type";
    }
    return $class;
}

sub report_stats {
    my $curr_chr = shift;
    my $curr_pos = shift;
    my $lib_id2  = shift;
    my $curr_chunk   = shift;
    my $chunks   = shift;
    die "Somethiong is wrong with BAM file - no pairs found" unless $stats{'total'};
    die "Somethiong is wrong with BAM file" unless $stats{'nonclonal'};
    
    warn "--------- Progress report -----------------------\n";
    warn "Currently at ",$lib_id2," lib, ",$cfg{$lib_id2}{'name'},"\n";
    warn "Currently at chunk ",($curr_chunk+1)," / $chunks\n";
    warn "Currently at chr $curr_chr : $curr_pos\n";
    warn "BAM lines reads  : ",$lines_read*2,"\n";
    warn "Pairs processed  : $stats{'total'}\n";
    warn "Non-clonal pairs : $stats{'nonclonal'}\t",100*$stats{'nonclonal'}/$stats{'total'},"%\n" if $cfg{'Distribution'}{'remove_clonal'};
    warn "Consistent pairs : $stats{'consistent'}\t",100*$stats{'consistent'}/$stats{'nonclonal'},"%\n";
    warn "Balance ",$stats{'balance'}, "\n";
    my $inv = $stats{'HH'}+$stats{'hh'}+$stats{'TT'}+$stats{'tt'};
    warn "Inverted pairs   : $inv\t",100*$inv/$stats{'nonclonal'},"%\n";
    my $anti = $stats{'HT'}+$stats{'ht'};
    warn "Everted pairs    : $anti\t",100*$anti/$stats{'nonclonal'},"%\n";
    warn "Remote pairs     : $stats{'remote'}\t",100*$stats{'remote'}/$stats{'nonclonal'},"%\n";
}

sub save_distr {
    my $lib_no = shift;
    my $cumul = 0;
    foreach my $size (sort {$a<=>$b} keys %distr ) {
        $cumul += $distr{$size};
        print FDIST join ("\t", $size, $distr{$size}/$stats{'consistent'}, $cumul/$stats{'consistent'}), "\n";
    }
}

sub save_stats {
    my $lib_no = shift;
    my $inverted = $stats{'HH'}+$stats{'hh'}+$stats{'TT'}+$stats{'tt'};
    my $everted  = $stats{'HT'}+ $stats{'ht'};
    print FSTAT join( "\n", 'Library_no'."\t".$lib_no,
			    'Library_name'."\t".$cfg{$lib_no}{'name'},
			    'Library_type'."\t".$cfg{$lib_no}{'type'},
			    'Pairs_total'."\t". $stats{'total'},
			    'Pairs_nonclonal'."\t".$stats{'nonclonal'},
			    'Mapped_consistently'."\t".$stats{'consistent'},
			    'Mapped_consistently_perc'."\t".(100*$stats{'consistent'}/$stats{'nonclonal'}),
			    'Mapped_inverted'."\t".$inverted,
			    'Mapped_inverted_perc'."\t".(100*$inverted/$stats{'nonclonal'}),
			    'Mapped_everted'."\t".$everted,
			    'Mapped_everted_perc'."\t".(100*$everted/$stats{'nonclonal'}),
			    'Mapped_remote'."\t".$stats{'remote'},
			    'Mapped_remote_perc'."\t".(100*$stats{'remote'}/$stats{'nonclonal'}),
		    ),"\n";
}
