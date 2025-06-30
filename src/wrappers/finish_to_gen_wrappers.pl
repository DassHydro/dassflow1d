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
# Replace function calc_cost_and_gradients in dassflow1d/__init__.py
open(my $in, "<dassflow1d/__init__.py"); # Ouvrir le fichier "file" en mode "lecture"

my @copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

open(my $out, ">dassflow1d/__init__.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

my $step = 0;
my $line = 0;

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    $line = $line + 1;
    if ( $line == 1 )
    {
        if ( m/from __future__(.*)/i )
        {
            print $out "$_";
            print $out "import numpy as np\n";
            next;
            
        } else {
            print $out "import numpy as np\n";
        }
    }
    if ( $step == 0)
    {
        if ( m/^(\s*)def calc_cost_and_gradients(.*)grad\, (.*)/i )
        {
            print $out "$1def calc_cost_and_gradients$2$3\n";
            $step = 1;
            next;
        }
        else
        {
          print $out "$_";
          next;
        }
    }
    elsif ( $step == 1 )
    {
        if ( m/^(\s*)cost = calc_cost_and_gradients(.*)grad\, (.*)/i )
        {
            print $out "$1cost, grad = calc_cost_and_gradients$2$3\n";
            next;
        }
        if ( m/^(\s*)grad \: float array.*/i )
        {
            next;
        }
        if ( m/^(\s*)cost \: float(.*)/i )
        {
            print $out "$_";
            print $out "$1grad \: float array$2\n";
            next;
        }
        if ( m/^(\s*)cost \= _dassflow1d(.*)/i )
        {
            print $out "$1grad = np.zeros(ctrl.x.size)\n$_";
            next;
        }
        if ( m/^(\s*)return cost(.*)/i )
        {
            print $out "$1return cost, grad$2\n";
            $step = 0;
            next;
        }
        else
        {
          print $out "$_";
          next;
        }
    }
}
        
# Add function calc_cost_and_gradients_scipy_minimize
print $out "\ndef calc_cost_and_gradients_scipy_minimize(x, mdl, ctrl, obs, steady=False, logger=None):\n";
print $out "\n    \"\"\"\"\n";
print $out "    cost, grad = calc_cost_and_gradients_scipy_minimize(x, mdl, ctrl, obs)\n";
print $out "\n\n";    
print $out "    Parameters\n";
print $out "    ----------\n";
print $out "    x : Current control vector values\n";
print $out "    mdl : Model\n";
print $out "    ctrl : Control\n";
print $out "    obs : Observations\n";
print $out "\n    Returns\n";
print $out "    -------\n";
print $out "    cost : float\n";
print $out "    grad : float array\n";
print $out "\n";
print $out "    \"\"\"\n\n";
print $out "    if logger is None:\n";
print $out "        logger = print\n";
print $out "\n";
print $out "    ctrl.x[:] = x[:]\n";
print $out "    grad = np.zeros(x.size)\n";
print $out "    if steady:\n";
print $out "        model.te = model.ts\n";
print $out "    cost = _dassflow1d.f90wrap_calc_cost_and_gradients(grad=grad, mdl=mdl._handle, \\\n";
print $out "        ctrl=ctrl._handle, obs=obs._handle)\n";
print $out "\n";
print $out "    if mdl.status == 10:\n";
print $out "        logger(\"Simulation failed: NaN value(s) detected at downstream boundary during standard step\")\n";
print $out "    elif mdl.status == 11:\n";
print $out "        logger(\"Simulation failed: NaN value(s) detected inside a segment during standard step\")\n";
print $out "    elif mdl.status == 20:\n";
print $out "        logger(\"Simulation failed: NaN value(s) detected on hydraulic variables during preissman timestep\")\n";
print $out "    elif mdl.status == 21:\n";
print $out "        logger(\"Simulation failed: NaN value(s) for h or Q during preissman timestep\")\n";
print $out "    elif mdl.status == 22:\n";
print $out "        logger(\"Simulation failed: Negative Debord operand detected  during preissman timestep\")\n";
print $out "    elif mdl.status != 0:\n";
print $out "        logger(\"Simulation failed: Simulation failed (status=\%i)\" \% mdl.status)\n";
print $out "\n";
print $out "    if mdl.status == 0:\n";
print $out "        if mdl.warning_counters[0] > 0:\n";
print $out "            # heps correction if counted twice (direct and reverse mode)\n";
print $out "            logger(\"heps correction applied \%i times\" \% (mdl.warning_counters[0] / 2))\n";
print $out "\n";
print $out "    mdl.__last_cost__ = cost\n";
print $out "    mdl.__last_grad__ = grad\n";
print $out "    return cost, grad\n";
print $out "\n";
        
# Add function calc_cost_gradients_and_status
print $out "\ndef calc_cost_gradients_and_status(x, mdl, ctrl, obs, steady=False):\n";
print $out "\n    \"\"\"\"\n";
print $out "    cost, grad, status = calc_cost_gradients_and_status(x, mdl, ctrl, obs)\n";
print $out "\n\n";    
print $out "    Parameters\n";
print $out "    ----------\n";
print $out "    x : Current control vector values\n";
print $out "    mdl : Model\n";
print $out "    ctrl : Control\n";
print $out "    obs : Observations\n";
print $out "\n    Returns\n";
print $out "    -------\n";
print $out "    cost : float\n";
print $out "    grad : float array\n";
print $out "\n";
print $out "    \"\"\"\n\n";
print $out "    ctrl.x[:] = x[:]\n";
print $out "    grad = np.zeros(x.size)\n";
print $out "    if steady:\n";
print $out "        model.te = model.ts\n";
print $out "    cost = _dassflow1d.f90wrap_calc_cost_and_gradients(grad=grad, mdl=mdl._handle, \\\n";
print $out "        ctrl=ctrl._handle, obs=obs._handle)\n";
print $out "    mdl.__last_cost__ = cost\n";
print $out "    mdl.__last_grad__ = grad\n";
print $out "    return cost, grad, mdl.status\n";
print $out "\n";

# Add version
print $out "\n__version__ = '2.1.0'\n";
print $out "\n";

# Close file
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add function get_item_slice in dassflow1d/m_control.py
open($in, "<dassflow1d/m_control.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_control.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Control\((.*)/i )
        {
            $step = 1;
        }
        print $out "$_";
        next;
    }
    else
    {
        if (m/^(\s*)def temporal_correlation_array_fixed\(self, t, mu, n, m, array\):(.*)/i)
        {
            print $out "    def temporal_correlation_array(self, t, mu):\n";
            print $out "        array = numpy.zeros((t.size, t.size), order='F')\n";
            print $out "        self.temporal_correlation_array_fixed(t, mu, t.size, t.size, array)\n";
            print $out "        return array\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        elsif (m/^(\s*)def get_spatial_correlation_array_size\(self, msh\):(.*)/i)
        {
            print $out "    def spatial_correlation_array(self, msh, mu):\n";
            print $out "        size_bn = self.get_spatial_correlation_array_size(msh)\n";
            print $out "        array = numpy.zeros((size_bn, size_bn), order='F')\n";
            print $out "        self.spatial_correlation_array_fixed(msh, mu, size_bn, size_bn, array)\n";
            print $out "        return array\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        elsif (m/^(\s*)def get_field_correlation_array_size\(self, msh, field\):(.*)/i)
        {
            print $out "    def field_correlation_array(self, msh, field, mu):\n";
            print $out "        size_bn = self.get_field_correlation_array_size(msh, field)\n";
            print $out "        array = numpy.zeros((size_bn, size_bn), order='F')\n";
            print $out "        self.field_correlation_array_fixed(msh, field, mu, size_bn, size_bn, array)\n";
            print $out "        return array\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        elsif (m/^(\s*)def __str__\(self\):(.*)/i)
        {
            print $out "    def get_item_slice(self, index):\n";
            print $out "        if index < 0 or index > self.get_items_count():\n";
            print $out "            raise IndexError('index out of range')\n";
            print $out "        return slice(self.items[index].offset, self.items[index].offset+self.items[index].nx)\n";
            print $out "\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add functions in dassflow1d/m_mesh.py
open($in, "<dassflow1d/m_mesh.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_mesh.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Mesh\((.*)/i )
        {
            $step = 1;
        }
        print $out "$_";
        if (m/^(\s*)import logging/i )
        {
            print $out "import numpy as np\n";
        }
        next;
    }
    else
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
    
            print $out "    def get_segment_field(self, iseg, field, base_cs=True):\n";
            print $out "    \n";
            print $out "        if field == \"curvilinear_abscissa\" or field == \"x\" or field == \"xs\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].x for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].x for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"coords\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([[self.cs[i].coord.x, self.cs[i].coord.y] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([[self.cs[i].coord.x, self.cs[i].coord.y] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"bathy\" or field == \"b\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].bathy for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].bathy for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"z0\" or field == \"zlow\" or field == \"Hlow\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].level_heights[0] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].level_heights[0] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"zm\" or field == \"zhigh\" or field == \"Hhigh\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].level_heights[-1] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].level_heights[-1] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"w0\" or field == \"wlow\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].level_widths[0] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].level_widths[0] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"zmax\" or field == \"Hmax\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].level_heights[self.cs[i].nlevels-1] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].level_heights[self.cs[i].nlevels-1] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"wmax\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].level_widths[self.cs[i].nlevels-1] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].level_widths[self.cs[i].nlevels-1] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"A0\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            w0 = [self.cs[i].level_widths[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "            z0 = [self.cs[i].level_heights[0] for i in range(first_cs-1, last_cs)]\n";
            print $out "            bathy = [self.cs[i].bathy for i in range(first_cs-1, last_cs)]\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([(z0[i]-bathy[i])*w0[i] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([(z0[i]-bathy[i])*w0[i] for i in range(0, len(w0))])\n";
            print $out "        elif field == \"kpar1\":\n";
            print $out "            first_cs = self.seg[iseg].first_cs\n";
            print $out "            last_cs = self.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                return np.array([self.cs[i].strickler_params[0] for i in indices])\n";
            print $out "            else:\n";
            print $out "                return np.array([self.cs[i].strickler_params[0] for i in range(first_cs-1, last_cs)])\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown segment field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add functions in dassflow1d/m_obs.py
open($in, "<dassflow1d/m_obs.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_obs.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        if (m/^(\s*)class Obsstation\((.*)/i )
        {
            $step = 1;
        }
        print $out "$_";
        if (m/^(\s*)import logging(.*)/i )
        {
            print $out "import numpy as np\n";
        }
        next;
    }
    elsif ($step == 1)
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
    
            print $out "    def setup(self, mesh, t, values, iseg=None, indices=None, coords=None, x=None):\n";
            print $out "    \n";
            print $out "        if indices is not None and coords is not None:\n";
            print $out "            raise ValueError(\"indices and coords cannot be specified at the same time\")\n";
            print $out "        if indices is not None and x is not None:\n";
            print $out "            raise ValueError(\"indices and curvilinear abscissa cannot be specified at the same time\")\n";
            print $out "        if coords is not None and x is not None:\n";
            print $out "            raise ValueError(\"coords and curvilinear abscissa cannot be specified at the same time\")\n";
            print $out "        if indices is not None:\n";
            print $out "            self.setup_station_using_indices(mesh, indices, t, values)\n";
            print $out "        elif coords is not None:\n";
            print $out "            if isinstance(coords, tuple):\n";
            print $out "                coords = np.array(list(coords)).reshape((2, 1))\n";
            print $out "            elif isinstance(coords, list):\n";
            print $out "                coords = np.array(coords)\n";
            print $out "                if coords.ndim == 1:\n";
            print $out "                    coords = coords.reshape((2, 1))\n";
            print $out "                elif coords.shape[2] == 2:\n";
            print $out "                    coords = coords.T\n";
            print $out "            if iseg is not None:\n";
            print $out "                self.setup_station_using_coords_and_segment(mesh, iseg, coords, t, values)\n";
            print $out "            else:\n";
            print $out "                self.setup_station_using_coords(mesh, coords, t, values)\n";
            print $out "        elif x is not None:\n";
            print $out "            if iseg is not None:\n";
            print $out "                self.setup_station_using_abscissa_and_segment(mesh, iseg, x, t, values)\n";
            print $out "            else:\n";
            print $out "                self.setup_station_using_abscissa(mesh, x, t, values)\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"indices or coords or curvilinear abscissa must be specified\")\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;

#-----------------------------------------------------------------------------------------------------------------------
# Add functions in dassflow1d/m_sw_mono.py
open($in, "<dassflow1d/m_sw_mono.py"); # Ouvrir le fichier "file" en mode "lecture"

@copy = <$in>;  #Copier le contenu du fichier dans "copy"

close $in; #Fermer l'original

$step = 0;

open($out, ">dassflow1d/m_sw_mono.py"); # Reouvrir le fichier en mode "ecriture"  (écrasement du fichier)

for ( @copy ) #Parcourt la copie du fichier ligne par ligne
{
    if ( $step == 0)
    {
        
        if (m/^(\s*)class Unknowns\((.*)/i )
        {
            $step = 1;
        }
        if (m/^(\s*)class Model\((.*)/i )
        {
            $step = 2;
        }
        print $out "$_";
        if (m/^(\s*)import logging/i )
        {
            print $out "import numpy as np\n";
        }
        next;
    }
    elsif ($step == 1)
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
    
            print $out "    def get_segment_field(self, mesh, iseg, field, base_cs=True):\n";
            print $out "    \n";
            print $out "        if field == \"coords\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(mesh.cs[i].coords)\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([mesh.cs[i].coords for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"a\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.a[i])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.a[i] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"q\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.q[i])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.q[i] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"h\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.h[i])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.h[i] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"z\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.h[i] + mesh.cs[i].bathy)\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.h[i] + mesh.cs[i].bathy for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"zlow\":\n";
            print $out "            first_cs = mesh.seg[iseg].first_cs\n";
            print $out "            last_cs = mesh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if mesh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(mesh.cs[i].level_heights[0])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([mesh.cs[i].level_heights[0] for i in range(first_cs-1, last_cs)])\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown segment field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
    elsif ($step == 2)
    {
        if (m/^(\s*)def __str__\(self\):(.*)/i)
        {
            print $out "    def get_results_time(self, iout):\n";
            print $out "    \n";
            print $out "        return self.res.t[iout]\n";
            print $out "\n";
            print $out "    def get_segment_results(self, iseg, iout, field, base_cs=True):\n";
            print $out "    \n";
            print $out "        if field == \"q\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.res.q[i, iout])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.res.q[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"h\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.res.h[i, iout])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.res.h[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"z\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.res.h[i, iout] + self.msh.cs[i].bathy)\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.res.h[i, iout] + self.msh.cs[i].bathy for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"w\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.msh.cs[i].htow_noupdate(self.res.h[i, iout]))\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                w = np.array([self.msh.cs[i].htow_noupdate(self.res.h[i, iout]) for i in range(first_cs-1, last_cs)])\n";
            print $out "            return w\n";
            print $out "        elif field == \"a\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.res.a[i, iout])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.res.a[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"p\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(self.msh.cs[i].htop(self.res.h[i, iout]))\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return np.array([self.msh.cs[i].htop(self.res.h[i, iout]) for i in range(first_cs-1, last_cs)])\n";
            print $out "        elif field == \"Fr\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            if base_cs is True:\n";
            print $out "                indices = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        indices.append(i)\n";
            print $out "                w = np.array([self.msh.cs[i].htow_noupdate(self.res.h[i, iout]) for i in indices])\n";
            print $out "                q = np.array([self.res.q[i, iout] for i in indices])\n";
            print $out "                a = np.array([self.res.a[i, iout] for i in indices])\n";
            print $out "            else:\n";
            print $out "                w = np.array([self.msh.cs[i].htow_noupdate(self.res.h[i, iout]) for i in range(first_cs-1, last_cs)])\n";
            print $out "                q = np.array([self.res.q[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "                a = np.array([self.res.a[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "            return np.sqrt(q**2 * w / (self.gravity * a**3))\n";
            print $out "        elif field == \"s\":\n";
            print $out "            first_cs = self.msh.seg[iseg].first_cs\n";
            print $out "            last_cs = self.msh.seg[iseg].last_cs\n";
            print $out "            h = np.array([self.res.h[i, iout] for i in range(first_cs-1, last_cs)])\n";
            print $out "            s = np.zeros(last_cs - first_cs + 1)\n";
            print $out "            self.free_surface_slopes_segment(iseg+1, h, s)\n";
            print $out "            if base_cs is True:\n";
            print $out "                values_list = []\n";
            print $out "                for i in range(first_cs-1, last_cs):\n";
            print $out "                    if self.msh.cs[i].ibase > 0:\n";
            print $out "                        values_list.append(s[i - first_cs + 1])\n";
            print $out "                return np.array(values_list)\n";
            print $out "            else:\n";
            print $out "                return s\n";
            print $out "        else:\n";
            print $out "          raise ValueError(\"Unknown result field: \%s\" \% field)\n";
            print $out "\n";
            print $out "$_";
            next;
        }
        else
        {
          if (m/^(\s*)class (.*)/i)
          {
              $step = 0;
          }
          print $out "$_";
          next;
        }
    }
}
close $out;


#-----------------------------------------------------------------------------------------------------------------------
open($out, ">env.sh"); # ouvrir le fichier en mode "ecriture"  (écrasement du fichier)

print $out "FILEPATH=`realpath \${BASH_SOURCE[0]}`\n";
print $out "DIR=`dirname \$FILEPATH`\n";
print $out "export PYTHONPATH=\$DIR:\$PYTHONPATH\n";

close($out);


#-----------------------------------------------------------------------------------------------------------------------
open($out, ">../env.sh"); # ouvrir le fichier en mode "ecriture"  (écrasement du fichier)

print $out "FILEPATH=`realpath \${BASH_SOURCE[0]}`\n";
print $out "DIR=`dirname \$FILEPATH`\n";
print $out "export DASSFLOW1D_BINDIR=\$DIR/bin\n";
print $out "export PYTHONPATH=\$DIR/bin:\$DIR/api:\$PYTHONPATH\n";
print $out "export PATH=\$DIR/bin:\$PATH\n";

close($out);

