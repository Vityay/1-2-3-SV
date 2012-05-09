#!/usr/bin/perl -w
use strict;
use Config::Std;

# 1-2-3-SV step 1
# v 0.04
# Last modified 20120506

STDOUT->autoflush(1);

my %cfg = ();
my ( $chunk, $lib_id );
my $chunk_count = 1;
my $f = 0;
my $r = 0;
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
my %nuq = ();

foreach $lib_id ( sort {$a<=>$b} keys %lib_ids ) { 

    if (-e 'sv_size_distribution.'.$cfg{$lib_id}{'name'} ) {
         warn "lib $lib_id seems to be processed, skipping it\n";
         next;
    }

    die 'Incorrect number of chunks is given for library ',$lib_id, "\n" if $cfg{$lib_id}{'chunks'} !~ m/^[1-9]\d*$/;

    %distr = ();
    %stats = ();
    %seen = ();
    %fwd = ();
    %rev = ();
    $chunk_count = $cfg{$lib_id}{'chunks'};
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
    $lines_read = 0;

    my $lib = $cfg{$lib_id}{'name'};
    my $libprefix = $cfg{$lib_id}{'prefix'};
    my $lib_type  = uc( $cfg{$lib_id}{'type'} ) or die "Library type for library $lib_id not specified";
    warn 'Determining insert sizes for lib # ', $lib_id, ' : type ', $lib_type, ' in ', $cfg{$lib_id}{'chunks'}, " chunks\n";
    open FDIST, '>', 'sv_size_distribution.'.$lib;
    open FSTAT, '>', 'sv_pairs_stats.'.$lib;
    open FNUCR, '>', 'sv_pairs_remote.'.$lib;
    open FNUCI, '>', 'sv_pairs_inversions.'.$lib;
    open FNUCE, '>', 'sv_pairs_evertions.'.$lib;
    foreach my $chunk2 ( 0 .. ( $cfg{$lib_id}{'chunks'} - 1 ) ) {
        $chunk = $chunk2;
        undef %fwd;
        undef %rev;
        warn 'Started with chunk ', ($chunk + 1), ' of ', $cfg{$lib_id}{'chunks'}, ' for lib # ', $lib_id,' : ', $lib,' type ', $lib_type, "\n";
        open F1, $cfg{'General'}{'samtools_exe'}.' view -F 4 '.$cfg{$lib_id}{'fileF'}.'|' or die 'Could not open file with forward tags for lib '.$lib;
        open F2, $cfg{'General'}{'samtools_exe'}.' view -F 4 '.$cfg{$lib_id}{'fileR'}.'|' or die 'Could not open file with reverse tags for lib '.$lib;
        while ( !eof(F1) or !eof(F2) ) {
            # Forward read
            unless (eof(F1)) {
                my $line1 = <F1>;
                my ( $clone, $chr, $pos, $strand, $dir ) = get_coord( $lib_id, $line1 );
                report_stats( $chr, $pos, $lib_id, $chunk, $cfg{$lib_id}{'chunks'} ) if $chr and not(++$lines_read % 10_000);
                if ( !$clone or not(in_chunk($clone)) ) { # not from lib or chunk
                    ; #do nothings
                }
                elsif ( $clone and !defined($chr) ) { # unmapped or non-unique
		    if ( exists($fwd{$clone}) or exists($rev{$clone}) or exists($nuq{$clone}) ) { # erase if pair seen
		        free($clone);
		    }
		    else { #mark if pair not seen
		         $nuq{$clone} = 1;
		    }
                }
                else { # mapped uniquely
                    if ( exists($nuq{$clone}) ) {
                        free($clone);
                    }
                    else {
                        if ( $dir eq 'F' ) {
                            $fwd{$clone} = join("\t", $chr, $strand, $pos );
                            $f++; 
                            process_ditag( $clone, $lib_type ) if ( exists( $rev{$clone} ) );
                        }
                        elsif ($dir eq 'R' ) {
                            $rev{$clone} = join("\t", $chr, $strand, $pos );
                            $r++; 
                            process_ditag( $clone, $lib_type ) if ( exists( $fwd{$clone} ) );
                        }
                    }
                }
            }

            # Reverse read
            unless ( $cfg{$lib_id}{'fileR'}and eof(F2) ) {
                my $line1 = <F2>;
                my ( $clone, $chr, $pos, $strand, $dir ) = get_coord( $lib_id, $line1 );
                $dir = 'R';
                #report_stats( $chr, $pos, $lib_id, $chunk, $cfg{$lib_id}{'chunks'} ) if $chr and not( ++$lines_read % 10_000);
                if ( !$clone or not(in_chunk($clone)) ) { # not from lib or chunk
                    ; #do nothing
                }
                elsif ( $clone and !defined($chr) ) { # unmapped or non-unique
		    if ( exists($fwd{$clone}) or exists($rev{$clone}) or exists($nuq{$clone}) ) { # erase if pair seen
		        free($clone);
		    }
		    else { #mark if pair not seen
		         $nuq{$clone} = 1;
		    }
                }
                else { # mapped uniquely
                    if ( exists($nuq{$clone}) ) {
                        free($clone);
                    }
                    else {
                        $rev{$clone} = join("\t", $chr, $strand, $pos );
                        $r++; 
                        process_ditag( $clone, $lib_type ) if ( exists( $fwd{$clone} ) );
                    }
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
         return() unless $clone =~ m/$prefix/;
    }
    my $trim = $cfg{$lib}{'trim'} if exists($cfg{$lib}{'trim'});
    $clone =~ s/$trim// if defined($trim);

    if ( $map_flag & 4 ) {
        free( $clone );
        return($clone);
    }

    return($clone) if exists($cfg{'Distribution'}{'max_ambiguity'} ) and $line =~ m/X0\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'max_ambiguity' };
    return($clone) if exists($cfg{'Distribution'}{'secondary_hits'}) and $line =~ m/X1\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'secondary_hits'};
    return($clone) if exists($cfg{'Distribution'}{'alignment_gaps'}) and $line =~ m/XO\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'alignment_gaps'};
    return($clone) if exists($cfg{'Distribution'}{'alignment_mismatches'}) and $line =~ m/XM\:i\:(\d+)/ and $1 > $cfg{'Distribution'}{'alignment_mismatches'};

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
    my $dir = $map_flag & 128 ? 'R' : 'F';
    return( $clone, $c, $p, $strand, $dir );
}

sub process_ditag {
    my $clo = shift;
    my $lib_tp = shift;
    my ( $c1, $s1, $p1 ) = split /\t/, $fwd{$clo};
    my ( $c2, $s2, $p2 ) = split /\t/, $rev{$clo};

    my $ab = 1;
    $stats{'total'}++;
    if ( $c1 ne $c2 or abs($p1 - $p2 ) > $cfg{'Distribution'}{'max_distance'} ) { # remotely located di-tags
        return if $cfg{'Distribution'}{'remove_clonal'} and ++$seen{ join( "\t", $c1, $s1, $p1, $c2, $s2, $p2 )} > 1; # we do not take it into account 
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
            return if $cfg{'Distribution'}{'remove_clonal'} and ++$seen{ join( "\t", $c1, $s1, $p1, $c2, $s2, $p2 )} > 1; # we do not take it into account 
            ( $p1, $p2 ) = ( $p2, $p1 ) if $p1 > $p2; 
            print FNUCE join( "\t", $c1, int($p1/$cfg{'General'}{'bin_size'}), int($p2/$cfg{'General'}{'bin_size'})), "\n";
        }
        else { #inversion
            return if $cfg{'Distribution'}{'remove_clonal'} and ++$seen{ join( "\t", $c1, $s1, $p1, $c2, $s2, $p2 )} > 1; # we do not take it into account 
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
    my ( $curr_chr, $curr_pos, $lib_id2, $curr_chunk, $chunks ) = @_;
    die "Somethiong is wrong with BAM file - no pairs found" unless $stats{'total'};
    die "Somethiong is wrong with BAM file" unless $stats{'nonclonal'};
    my $inv = 0;
    foreach ( 'HH', 'hh', 'TT', 'tt' ) {
        $inv += $stats{$_} if exists( $stats{$_} );
    }    
    my $evr = 0;
    foreach ( 'HT', 'ht' ) {
        $evr += $stats{$_} if exists( $stats{$_} );
    }
    my $rem = 0;
    $rem+=$stats{'remote'} if exists( $stats{'remote'} );
    print join( ' ', 'now at', $curr_chr.':'.$curr_pos,
                    'pairs:'.$stats{'total'},
                    'good:'.$stats{'consistent'},
                    'inv:'.$inv,
                    'evr:'.$evr,
                    'rem:'.$rem,
                    '       ',
              ), "\r";
    clean($lib_id2) if $lines_read % 10_000_000 == 0 and scalar(keys %fwd) + scalar(keys %rev) + scalar(keys %nuq) > 10_000_000;
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
			    'keys fwd'."\t".scalar(keys %fwd),
			    'keys rev'."\t".scalar(keys %rev),
			    'keys nuq'."\t".scalar(keys %nuq),
			    'keys seen'."\t".scalar(keys %seen),
		    ),"\n";
}

sub free {
    my $clone = shift;
    delete $fwd{$clone} if exists($fwd{$clone});
    delete $rev{$clone} if exists($rev{$clone});
    delete $nuq{$clone} if exists($nuq{$clone});
}

sub in_chunk {
    my $str = shift;
    my ( $sum, $cnt ) = ( 0, 0 );
    map { $sum += ord($_)*(++$cnt) } split //, $str;
    return 1 if $sum % $chunk_count == $chunk;
    return 0;
}

sub clean {
    my $lib_no = shift;
    my $trim = $cfg{$lib_no}{'trim'} if exists($cfg{$lib_no}{'trim'});
    my @files = ( $cfg{$lib_no}{'fileF'} );
    push @files, $cfg{$lib_no}{'fileR'} if $cfg{$lib_no}{'fileR'};
    print "Cleaning buffers\n";
    foreach my $file ( @files ) {
        my $i = 0;
        open F3, $cfg{'General'}{'samtools_exe'}.' view -f 4 '.$file.' |' or die 'Could not open file '.$file.' for lib '.$lib_no;
        while ( <F3> ) {
            my ( $clone, $map_flag, $c, $p, $qual, $match ) = split /\t/;
  
            # Skip if prefix unless multi-lib BAM and prefix is specified
            if ( $cfg{$lib_no}{'prefix'} ) {
                 my $prefix = $cfg{$lib_no}{'prefix'};
                 next unless $clone =~ m/$prefix/;
            }
            $clone =~ s/$trim// if defined($trim);
            delete $fwd{$clone} if exists($fwd{$clone});
            delete $rev{$clone} if exists($rev{$clone});
            delete $nuq{$clone} if exists($nuq{$clone});
            $i++;
        }
        close F3;
    }
}