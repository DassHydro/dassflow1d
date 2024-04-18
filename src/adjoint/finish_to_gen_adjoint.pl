#======================================================================================================================!
#
#                    DassFlow Version 2.0
#
#======================================================================================================================!
#
#  Copyright University of Toulouse-INSA & CNRS (France)
#
#  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
#  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
#  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
#  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
#
#  DassFlow software includes few mostly independent "modules" with common architectures and structures:
#    - Shallow Module (Shallow Water Model, Finite Volume Method), i.e. the present code.
#    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
#  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
#
#  Many people have contributed to the DassFlow development from the initial version to the latest ones.
# 	Current main developer:
#               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
# 	with scientific and/or programming contributions of:
#               R. Madec   (Mathematics Institute of Toulouse IMT).
#               K. Larnier (Fluid Mechanics Institute of Toulouse IMFT).
#               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
#               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
#	and former other developers (M. Honnorat and J. Marin).
#
#  Scientific Contact : jerome.monnier@insa-toulouse.fr
#  Technical  Contact : frederic.couderc@math.univ-toulouse.fr
#
#  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
#  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
#  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
#  license, users are provided only with a limited warranty and the software's author, the holder of the economic
#  rights, and the successive licensors have only limited liability.
#
#  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
#  developing or reproducing the software by the user in light of its specific status of free software, that may
#  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
#  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
#  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
#  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
#  accept its terms.
#
#======================================================================================================================!

#!/usr/bin/env perl

use strict;
use warnings;

#-----------------------------------------------------------------------------------------------------------------------
# Lecture de l'intégralité des fichiers du répertoire ./tap qui sont du type 'm_*_diff.f90' et 'm_*_back.f90'
#  - Le fichier est conservé si et seulement si il contient au moins une subroutine ou fonction différenciée
#  - Tous les fichiers ayant une subroutine différenciée sont agrégés à la liste  list
#  - Supprime les subroutines et les fonctions non différenciées
#-----------------------------------------------------------------------------------------------------------------------

my $debug = 0;
my @list; #Initialisation de la variable list

foreach my $file ( <m_*_diff.f90> , <m_*_back.f90> ) # Parmi tous les fichiers de type 'm_*_diff.f90' et 'm_*_back.f90'
   {
   
   open(my $in, "<$file"); # Ouvrir le fichier "file" en mode "lecture"

   my @copy = <$in>;  #Copier le contenu du fichier dans "copy"

   close $in; #Fermer l'original

   open(my $out, ">$file"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

   my $test = 0; # Initialisation variable "test" à 0
   my $in_back_type = 0; # Initialisation variable "test" à 0
   my $keep = 0; # Initialisation variable "keep" à 0
   
   if ( $debug == 1 )
   {
      print "----- FILE IS $file -----\n";
   }
 
   for ( @copy ) #Parcourt la copie du fichier ligne par ligne
   {
        
      if ( $debug == 1 )
      {
         print "(LINE)$_\n(test=$test)\n";
      }
      
      if ($test == 2)
      {
         if( m/^\s*END\s*SUBROUTINE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
         {
            print $out $_; # On écrit la ligne dans le fichier out
            $test = 0;
            next;
         }
         else
         {
            print $out $_; # On écrit la ligne dans le fichier out
            next;
         }
      }
      if ( m/^(\s*)TYPE\((\w*)_(DIFF|BACK)\)(.*)/i ) 
         {
         # Nothing todo
         }

      elsif    ( m/^\s*TYPE (\w*)_(DIFF|BACK).*/i ) #Si la copie contient un type differencie
      {
         $in_back_type = 1;
         $keep = 1;     #On fixe keep à 1
      }
      elsif ( m/^\s*TYPE.*/i )  #Sinon si la copie contient une surboutine
      {
            if ($in_back_type == 0)
            {
                $test = 1;     #On fixe test à 1
                print $out ""; 
            }
            else
            {
                print $out $_; 
            }
         next;
      }
      if ( m/^(\s*)USE (.*)_(DIFF|BACK), ONLY(.*),/i )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
      {
           
         my $sub1 = lc $3;
         print $out "$1USE $2, ONLY $4 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine
         print $out "$1USE $2_$3, ONLY $4_$sub1 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine

#          my $sub1 = $1;
#          my $sub2 = $2;
#          my $sub3 = $3;
# 
#          if ( "@list_sort" =~ /.*$sub2.*/i ) 
#             {
#             print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
#             }

         next;
      }
      if    ( m/^\s*SUBROUTINE (\w*)_(DIFF|BACK).*/i ) #Si la copie contient une subroutine differenciee 
      {
         $test = 2;
         if ( $1 !~ m/com_dof/i ) #Si la subroutine n'est pas com_dof
         {
            $keep = 1;              #On fixe keep à 1
            @list = (@list,$file);  #On ajoute le fichier à la liste
         }
      }

      elsif ( m/^\s*SUBROUTINE.*/i )  #Sinon si la copie contient une surboutine
      {
         $test = 1;     #On fixe test à 1
         print $out ""; 
         next;
      }
      
      if    (( m/^\s*.*FUNCTION\s*.*_(DIFF|BACK).*/i )&&( m/^\s.[^'END']\s*.*/i))
      {
         $keep = 1;     #On fixe keep à 1
      }

      elsif (( m/^\s*.*FUNCTION\s*.*/i )&&( m/^\s.[^'END']\s*.*/i))
      {
         $test = 1;     # On fixe test à 1
         print $out ""; 
         next;
      }

      if    ( m/^\s*END\s*TYPE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
      {
         $in_back_type = 0;
         print $out $_; # On écrit la ligne dans le fichier out
         next;
      }

      elsif ( m/^\s*END\s*TYPE.*/i ) #Si la copie contient une fin de subroutine 
      {
         $test = 0;     #On fixe test à 0
         print $out ""; 
         next;
      }

      if    ( m/^\s*END\s*SUBROUTINE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
      {
         print $out $_; # On écrit la ligne dans le fichier out
         next;
      }

      elsif ( m/^\s*END\s*SUBROUTINE.*/i ) #Si la copie contient une fin de subroutine 
      {
         $test = 0;     #On fixe test à 0
         print $out ""; 
         next;
      }

      if    ( m/^\s*END\s*FUNCTION\s*.*_(DIFF|BACK).*/i ) #Si la copie  contient une fin de fonction differenciee
      { 
         print $out $_; # On écrit la ligne dans le fichier out
         next;
      }

      elsif ( m/^\s*END\s*FUNCTION\s*.*/i ) # Si la copie contient une fin de fonction
      {
         $test = 0;        # On fixe test à 0
         print $out "";    
         next;
      }

      if ( m/^\s*MODULE\s*([a-zA-Z0-9_]*)_(DIFF|BACK).*/i ) # Si la copie contient un module différencié
      {
         $test = 1;                                        # On fixe le test à 1
#          print $out $_."\n  USE $1\n  USE M_TAP_VARS\n\n"; # On écrit la ligne de déf du module plus l'inclusion de mpi et tap_vars
         print $out $_."\n  USE $1\n"; # On écrit la ligne de déf du module plus l'inclusion de mpi et tap_vars
         #print $out $_."\n  USE $1\n  USE M_MPI_$2\n  USE M_TAP_VARS\n\n";
         next;
      }

      if    ( m/^\s*END\s*MODULE\s*.*_(DIFF|BACK).*/i ) #Si la copie  contient une fin de fonction differenciee
      { 
         print $out $_; # On écrit la ligne dans le fichier out
         next;
      }

      if ( m/^\s*CONTAINS.*/i ) # Si la copie contient un contains
      {
         $test = 0; # On fixe test ) 0
      }

      if ( $test == 1 ) # Si la variable test est à 1
      {
         print $out ""; # On n'est écrit rien
         next;
      }

      print $out $_;  #On écrit la ligne dans le fichier out

      }

   close $out; # On ferme le fichier out

   if ( $keep == 0 ) #Si la variable keep est égal à 0 
      {
      unlink $file   #On supprime le fichier
      }

   }

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

my (%saw,@list_sort)=();
undef %saw;
@list_sort = sort(grep(!$saw{$_}++, @list));

open(my $diffvars, ">diffvars.txt");   #Ouvre le fichier diffvars en mode écriture (écrasement des données)

foreach my $file ( <*_diff.f90> , <*_back.f90>,<*_back0.f90> )  # Parmi tous les fichiers de type '*_diff.f90' et '*_back.f90'
   {
   open(my $in, "<$file"); # Ouvrir le fichier "file" en mode "lecture"
   my @copy = <$in>;  #Copier le contenu du fichier dans "copy"
   close $in; #Fermer l'original
   open(my $out, ">$file"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

   for ( @copy ) #Parcourt la copie du fichier ligne par ligne
      {
      
      if ( m/^!.*Hint:\s*([\w%]*) .* \**([\w%]*)/i ) #Si l'occurence '!  Hint:' est trouvée
         {
         my $sub1 = $1; # Réccupération du tableau
         my $sub2 = $2; # Réccupération du nom du tableau

         if ( $sub1 !~ /ISIZE(.)OF.*/i ) #Si ce n'est pas un problème de taille de tableau
            {
            next; #On passe
            }

#          if ( $sub2 == "mdl\%msh\%strickler_fields" ) #Traitement special pour mdl%msh%strickler_fields
#             {
#             $sub2 = "mdl_back\%msh\%strickler_fields"; #On passe
#             }

         # Corrected following line (fomerly evaluates to 1 has it is a match condition result, 
         # not extration
         # my $sub3 = ( $sub1 =~ /ISIZE(.)OF.*/i );
         my ($sub3) = ( $sub1 =~ /ISIZE(.)OF.*/i );
         print $diffvars "$file\t\t$sub1\t\t$sub2\n";

         for ( @copy ) 
            {
            if ( !( m/^!.*Hint:\s*([\w%]*) .* \**([\w%]*)/i )) #Si l'occurence '!  Hint:' n'est pas trouvée
               {
               s/$sub1\b/size($sub2,$sub3)/g;
               }
            }

         next;
         }
      # Exemple !  Hint: ISIZE1OFDrfdof_v should be the size of dimension 1 of array *dof%v
      #$1 : ISIZE1OFDrfdof_v
      #$2 : *dof%v

#       if ( m/^(\s*)DO ii1=1,size(mdl%msh%strickler_fields,1).*/ && $file == "calc_cost_back.f90" )
#       {
#         print "FOUND!!!!\n";
#         print $out "$1DO ii1=1,size(mdl_back%msh%strickler_fields,1)\n";
#       }
      if ( m/^(\s*)DO ii1=1,size\(mdl%msh%strickler_fields,1\).*/i)
      {
        print "FOUND0000 !\n";
      }
      if ( m/^(\s*)DO ii1=1,size\(mdl%msh%strickler_fields,1\).*/i && $file eq "calc_cost_back.f90" )
      {
        print "FOUND!!!!\n";
        print $out "$1DO ii1=1,size(mdl_back%msh%strickler_fields,1)\n";
        next;
      }
      if ( m/^(\s*)DO ii2=1,size\(mdl%msh%seg%strickler_fields,1\).*/i && $file eq "calc_cost_back.f90" )
      {
        print "FOUND 01!!!!\n";
        print $out "$1DO ii2=1,size(mdl_back%msh%seg(ii1)%strickler_fields,1)\n";
        next;
      }
      if ( m/^(\s*)DO ii2=1,size\(msh%seg%strickler_fields,1\).*/i && $file eq "apply_strickler_fields_back.f90" )
      {
#         print "FOUND 02!!!!\n";
        print $out "$1DO ii2=1,size(msh%seg(ii1)%strickler_fields,1)\n";
        next;
      }


      if ( m/^\s*IMPLICIT NONE.*/i && $file !~ /m_.*/ ) # Recherche des IMPLICIT NONE dans les fichiers "m_..."
         {
#          print $out "\n  USE M_TAP_VARS ! Added by Perl Script -> Need to be filled !!!\n\n".$_; # Remplacement par USE M_TAP_VARS
         print $out $_; # Remplacement par USE M_TAP_VARS
         next;
         }

      if ( m/^(\s*)USE (.*)_(DIFF|BACK), ONLY(.*),(.*)/i && $file !~ /m_.*/ )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
#       if ( m/^(\s*)USE (.*)_(DIFF|BACK)(.*)/i && $file !~ /m_.*/ )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
         {
           
#          print "$file : $_ ";
         my $sub1 = lc $5;
         if ( $debug == 1 )
         {
            print "$file : $_ // $1;$2;$3;$4;$5!";
         }
         print $out "$1USE $2, ONLY $4 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine
         print $out "$1USE $2_$3, ONLY :$sub1 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine

#          my $sub1 = $1;
#          my $sub2 = $2;
#          my $sub3 = $3;
# 
#          if ( "@list_sort" =~ /.*$sub2.*/i ) 
#             {
#             print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
#             }

         next;

         }
      if ( m/^(\s*)USE (.*)_(DIFF|BACK), ONLY(.*)/i && $file !~ /m_.*/ )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
         {
         #print "$file : $_ ";
#          print $out "$1USE $2 ! Replaced by Perl Script\n"; # Ajout de USE NOM module d'origine
         print $out "$_";

#          my $sub1 = $1;
#          my $sub2 = $2;
#          my $sub3 = $3;
# 
#          if ( "@list_sort" =~ /.*$sub2.*/i ) 
#             {
#             print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
#             }
# 
         next;

         }

      if ( m/^(\s*)USE (.*)_(DIFF|BACK)/i && $file !~ /m_.*/ )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
         {
         #print "$file : $_ ";
#          print $out "$1USE $2 ! Replaced by Perl Script\n"; # Ajout de USE NOM module d'origine
         print $out "$1USE $2_$3 ! Replaced by Perl Script\n"; # Ajout de USE NOM module d'origine

         my $sub1 = $1;
         my $sub2 = $2;
         my $sub3 = $3;

         if ( "@list_sort" =~ /.*$sub2.*/i ) 
            {
            print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
            }

         next;

         }

      if ( m/^(\s*)TYPE\((\w*)_(DIFF|BACK)\), DIMENSION\(:\), POINTER(.*)/i ) 
         {
         print $out "$1TYPE($2_$3), DIMENSION(:), ALLOCATABLE$4 ! Replaced by Perl Script(3*)\n";
         next;
         }

      if ( m/^(\s*)TYPE\((\w*)_(DIFF|BACK)\)(.*)/i ) 
         {
         print $out "$1TYPE($2_$3)$4 ! Replaced by Perl Script(3+)\n";
         next;
         }

      if ( m/^!.*Hint.*/i ||
           m/.*DIFFSIZES.*/i ) 
         {
         next;
         }

#      if ( m/^(.*)CALL(.*)_C_(DIFF|BACK)(.*)/i ) 
#         {
#         print $out "$1CALL$2$4 ! Replaced by Perl Script\n";
#         next;
#         }

# 
#       if ( m/^(.*)_C_(DIFF|BACK)(.*)/i ) 
#          {
#          #print "$1$3";
#          print $out "$1$3 ! Replaced by Perl Script(4)\n";
#          next;
#          }

      if ( m/^(.*)SUBROUTINE (.*)_CB(.*)/i ) 
         {
         print $out "$_";
         next;
         }

      if ( m/^(.*)FUNCTION (.*)_CB(.*)/i ) 
         {
         print $out "$_";
         next;
         }

     if ( m/^(.*)CALL(.*)_CB(.*)/i ) 
        {
        print $out "$1CALL$2$3 ! Replaced by Perl Script(5)\n";
        next;
        }

      if ( m/^(.*)_CB(.*)/i ) 
         {
         print $out "$1$2 ! Replaced by Perl Script(6)\n";
         next;
         }

      if ( m/^.*POINTER.*/i ) 
         {
         next;
         }

      print $out $_;

   }

   close $out;

}

close $diffvars;

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

foreach my $file ( <*.f90> , <*.f90> ) {

   if ( $file =~ /.*mpi_(back|diff).*/ ) {
      unlink $file
   }

   if ( $file !~ /.*_(back|diff).*/ ) {
      unlink $file
   }

   if ( $file =~ /.*_c_(back|diff).*/ ) {
      unlink $file
   }

}

#-----------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------

foreach my $file ( <*_diff.f90> , <*_back.f90> ) {

   open(my $in, "<$file");

   my @copy = <$in>;

   close $in;

   for ( @copy ) {

      if( m/.*ISIZE.*OF.*/i ) {

         print "Missing in $file : $_";

      }

   }

}

#-----------------------------------------------------------------------------------------------------------------------
# Update time_loop_back.f90 (LSA, etc.)
#-----------------------------------------------------------------------------------------------------------------------
# open(my $in, "<time_loop_back.f90"); # Ouvrir le fichier "file" en mode "lecture"

# my @copy = <$in>;  #Copier le contenu du fichier dans "copy"

# close $in; #Fermer l'original

# open(my $out, ">time_loop_back.f90"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

# my $test = 0; # Initialisation variable "test" à 0
# my $in_back_type = 0; # Initialisation variable "test" à 0
# my $keep = 0; # Initialisation variable "keep" à 0

# if ( $debug == 1 )
# {
#    print "----- FILE IS time_loop_back.f90 -----\n";
# }

# for ( @copy ) #Parcourt la copie du fichier ligne par ligne
# {
      
#    if ( $debug == 1 )
#    {
#       print "(LINE)$_\n(test=$test)\n";
#    }
   
#    if ($test == 2)
#    {
#       if( m/^\s*END\s*SUBROUTINE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
#       {
#          print $out $_; # On écrit la ligne dans le fichier out
#          $test = 0;
#          next;
#       }
#       else
#       {
#          print $out $_; # On écrit la ligne dans le fichier out
#          next;
#       }
#    }
#    if ( m/^(\s*)TYPE\((\w*)_(DIFF|BACK)\)(.*)/i ) 
#       {
#       # Nothing todo
#       }

#    elsif    ( m/^\s*TYPE (\w*)_(DIFF|BACK).*/i ) #Si la copie contient un type differencie
#    {
#       $in_back_type = 1;
#       $keep = 1;     #On fixe keep à 1
#    }
#    elsif ( m/^\s*TYPE.*/i )  #Sinon si la copie contient une surboutine
#    {
#          if ($in_back_type == 0)
#          {
#                $test = 1;     #On fixe test à 1
#                print $out ""; 
#          }
#          else
#          {
#                print $out $_; 
#          }
#       next;
#    }
#    if ( m/^(\s*)USE (.*)_(DIFF|BACK), ONLY(.*),/i )  # Recherche des USE ..._DIFF ou ..._BACK dans les fichiers différents de "m_..."
#    {
         
#       my $sub1 = lc $3;
#       print $out "$1USE $2, ONLY $4 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine
#       print $out "$1USE $2_$3, ONLY $4_$sub1 ! Replaced by Perl Script(*)\n"; # Ajout de USE NOM module d'origine

# #          my $sub1 = $1;
# #          my $sub2 = $2;
# #          my $sub3 = $3;
# # 
# #          if ( "@list_sort" =~ /.*$sub2.*/i ) 
# #             {
# #             print $out $sub1 ; print $out "USE $sub2\_$sub3\n";
# #             }

#       next;
#    }
#    if    ( m/^\s*SUBROUTINE (\w*)_(DIFF|BACK).*/i ) #Si la copie contient une subroutine differenciee 
#    {
#       $test = 2;
#       if ( $1 !~ m/com_dof/i ) #Si la subroutine n'est pas com_dof
#       {
#          $keep = 1;              #On fixe keep à 1
#          @list = (@list,$file);  #On ajoute le fichier à la liste
#       }
#    }

#    elsif ( m/^\s*SUBROUTINE.*/i )  #Sinon si la copie contient une surboutine
#    {
#       $test = 1;     #On fixe test à 1
#       print $out ""; 
#       next;
#    }
   
#    if    (( m/^\s*.*FUNCTION\s*.*_(DIFF|BACK).*/i )&&( m/^\s.[^'END']\s*.*/i))
#    {
#       $keep = 1;     #On fixe keep à 1
#    }

#    elsif (( m/^\s*.*FUNCTION\s*.*/i )&&( m/^\s.[^'END']\s*.*/i))
#    {
#       $test = 1;     # On fixe test à 1
#       print $out ""; 
#       next;
#    }

#    if    ( m/^\s*END\s*TYPE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
#    {
#       $in_back_type = 0;
#       print $out $_; # On écrit la ligne dans le fichier out
#       next;
#    }

#    elsif ( m/^\s*END\s*TYPE.*/i ) #Si la copie contient une fin de subroutine 
#    {
#       $test = 0;     #On fixe test à 0
#       print $out ""; 
#       next;
#    }

#    if    ( m/^\s*END\s*SUBROUTINE.*_(DIFF|BACK).*/i ) #Si la copie contient une fin de subroutine differenciee
#    {
#       print $out $_; # On écrit la ligne dans le fichier out
#       next;
#    }

#    elsif ( m/^\s*END\s*SUBROUTINE.*/i ) #Si la copie contient une fin de subroutine 
#    {
#       $test = 0;     #On fixe test à 0
#       print $out ""; 
#       next;
#    }

#    if    ( m/^\s*END\s*FUNCTION\s*.*_(DIFF|BACK).*/i ) #Si la copie  contient une fin de fonction differenciee
#    { 
#       print $out $_; # On écrit la ligne dans le fichier out
#       next;
#    }

#    elsif ( m/^\s*END\s*FUNCTION\s*.*/i ) # Si la copie contient une fin de fonction
#    {
#       $test = 0;        # On fixe test à 0
#       print $out "";    
#       next;
#    }

#    if ( m/^\s*MODULE\s*([a-zA-Z0-9_]*)_(DIFF|BACK).*/i ) # Si la copie contient un module différencié
#    {
#       $test = 1;                                        # On fixe le test à 1
# #          print $out $_."\n  USE $1\n  USE M_TAP_VARS\n\n"; # On écrit la ligne de déf du module plus l'inclusion de mpi et tap_vars
#       print $out $_."\n  USE $1\n"; # On écrit la ligne de déf du module plus l'inclusion de mpi et tap_vars
#       #print $out $_."\n  USE $1\n  USE M_MPI_$2\n  USE M_TAP_VARS\n\n";
#       next;
#    }

#    if    ( m/^\s*END\s*MODULE\s*.*_(DIFF|BACK).*/i ) #Si la copie  contient une fin de fonction differenciee
#    { 
#       print $out $_; # On écrit la ligne dans le fichier out
#       next;
#    }

#    if ( m/^\s*CALL POPREAL8(mdl%tc).*/i ) # Si la copie contient un contains
#    {
#       print $out "#ifdef USE_LSA";
#       print $out "        if (lsa%enabled == .true.) then";
#       print $out "          if (int((mdl%tc - mdl%dt) / gsa%dtout) < int(mdl%tc / gsa%dtout)) then";
#       print $out "            call apply_control_back(gsa%ctrl, gsa%ctrl_back, mdl, mdl_back)";
#       print $out "            call write_lsa_results(mdl, gsa)";
#       print $out "          endif";
#       print $out "        endif";
#       print $out "#endif";
#       $applied1 = 1;
#    }

#    if ( $test == 1 ) # Si la variable test est à 1
#    {
#       print $out ""; # On n'est écrit rien
#       next;
#    }

#    if ( m/^\s*1,ad_count*/i )
#    {
#       print $out "#ifdef USE_LSA";
#       print $out "  if (lsa%enabled == .true.) then";
#       print $out "    call init_lsa_results(lsa)";
#       print $out "  endif";
#       print $out "#endif";
#    }

#    print $out $_;  #On écrit la ligne dans le fichier out

#    if ( m/^\s*ADMM_TAPENADE_INTERFACE.*/i )
#    {
#       print $out "#ifdef USE_LSA";
#       print $out "  USE LSA";
#       print $out "#endif";
#    }


# }

# close $out; # On ferme le fichier out
