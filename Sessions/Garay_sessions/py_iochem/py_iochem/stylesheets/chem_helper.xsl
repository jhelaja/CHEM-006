<?xml version="1.0" encoding="UTF-8"?>
<!--

    Create module - Create module inside the ioChem-BD software.
    Copyright © 2014 ioChem-BD (contact@iochem-bd.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:ckbk="http://my.safaribooksonline.com/book/xml/0596009747/numbers-and-math/77"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"
    xmlns:foldl-func="foldl-func"
    xmlns:f="http://fxsl.sf.net/"
    exclude-result-prefixes="xs f foldl-func"
    version="2.0">

    <xsl:variable name="AU_TO_EV" select="27.211" />
    <xsl:variable name="AU_TO_CM_1" select="219474" />
    <xsl:variable name="AU_TO_DEBYE" select="2.5415" />

    <xsl:variable name="EV_TO_AU" select="1 div $AU_TO_EV"/>
    <xsl:variable name="EV_TO_CM_1" select="8065.5446"/>
    <xsl:variable name="EV_TO_NM" select="1239.84187"/>

    <xsl:variable name="CM_1_TO_AU" select="1 div $AU_TO_CM_1"/>
    <xsl:variable name="CM_1_TO_EV" select="1 div $EV_TO_CM_1"/>

    <xsl:variable name="siUnitsPrefix" select="'^si$'"/>
    <xsl:variable name="nonsiUnitsPrefix" select="'^nonsi$'"/>
    <xsl:variable name="nonsi2UnitsPrefix" select="'^nonsi2$'"/>

    <xsl:variable name="siUnitsSymbolMap">
        <entry key="s">s</entry>
        <entry key="m">m</entry>
        <entry key="ampere">A</entry>
        <entry key="kg">kg</entry>
        <entry key="k">K</entry>
        <entry key="mol">mol</entry>
        <entry key="candela">cd</entry>
        <entry key="none">none</entry>
        <entry key="radian">rad</entry>
        <entry key="steradian">sr</entry>
        <entry key="hertz">Hz</entry>
        <entry key="newton">N</entry>
        <entry key="n.mol-1">N.mol-1</entry>
        <entry key="joule">J</entry>
        <entry key="watt">W</entry>
        <entry key="pascal">Pa</entry>
        <entry key="coulomb">C</entry>
        <entry key="volt">V</entry>
        <entry key="ohm">Ω</entry>
        <entry key="farad">F</entry>
        <entry key="siemens">S</entry>
        <entry key="weber">Wb</entry>
        <entry key="tesla">T</entry>
        <entry key="henry">H</entry>
        <entry key="becquerel">Bq</entry>
        <entry key="gray">Gy</entry>
        <entry key="sievert">Sv</entry>
        <entry key="katal">kat</entry>
        <entry key="mol.m-3.s-1">mol.m-3.s-1</entry>
        <entry key="mol.dm-3.s-1">mol.dm-3.s-1</entry>
        <entry key="m3.mol-1.s-1">m3.mol-1.s-1</entry>
        <entry key="m6.mol-2.s-1">m6.mol-2.s-1</entry>
        <entry key="molarity">_i_M__i_</entry>
        <entry key="molality">_i_m__i_</entry>
        <entry key="m2">m2</entry>
        <entry key="m3">m3</entry>
        <entry key="m-3">m-3</entry>
        <entry key="m-3.m-3">m-3.m-3</entry>
        <entry key="m-4">m-4</entry>
        <entry key="m.s-1">m.s-1</entry>
        <entry key="m.s-2">m.s-2</entry>
        <entry key="rad.s-1">rad.s-1</entry>
        <entry key="n.s">N.s</entry>
        <entry key="n.m.s">N.m.s</entry>
        <entry key="n.m">N.m</entry>
        <entry key="m-1">m-1</entry>
        <entry key="kg.m-3">kg.m-3</entry>
        <entry key="kg-1.m3">kg-1.m3</entry>
        <entry key="mol.mol-1">mol.mol-1</entry>
        <entry key="m-3.mol">m-3.mol</entry>
        <entry key="m3.mol-1">m3.mol-1</entry>
        <entry key="j.k-1">J.K-1</entry>
        <entry key="j.k-1.mol-1">J.K-1.mol-1</entry>
        <entry key="j.k-1.kg-1">J.K-1.kg-1</entry>
        <entry key="j.mol-1">J.mol-1</entry>
        <entry key="j.kg-1">J.kg-1</entry>
        <entry key="j.m-3">J.m-3</entry>
        <entry key="n.m-1">N.m-1 = J.m-2</entry>
        <entry key="w.m-2">W.m-2</entry>
        <entry key="w.m-1.k-1">W.m-1.K-1</entry>
        <entry key="m2.s-1">m2.s-1</entry>
        <entry key="m-2.s-1">m-2.s-1</entry>
        <entry key="mol.m-2.s-1">mol.m-2.s-1</entry>
        <entry key="kg.m-2.s-1">kg.m-2.s-1</entry>
        <entry key="pa.s">Pa.s = N.s.m-2</entry>
        <entry key="c.m-3">C.m-3</entry>
        <entry key="a.m-2">A.m-2</entry>
        <entry key="s.m-1">S.m-1</entry>
        <entry key="s.m2.mol-1">S.m2.mol-1</entry>
        <entry key="f.m-1">F.m-1</entry>
        <entry key="h.m-1">H.m-1</entry>
        <entry key="v.m-1">V.m-1</entry>
        <entry key="a.m-1">A.m-1</entry>
        <entry key="cd.m-2">cd.m-2</entry>
        <entry key="c.kg-1">C.kg-1</entry>
        <entry key="gy.s-1">Gy.s-1</entry>
        <entry key="j.m-1">J.m-1</entry>
    </xsl:variable>
    <xsl:variable name="nonsiUnitsSymbolMap">
        <entry key="degree">deg</entry>
        <entry key="degreePerMin">deg min-1</entry>
        <entry key="cm.s-1">cm.s-1</entry>
        <entry key="bohr">Bohr</entry>
        <entry key="angstrom">Å</entry>
        <entry key="dm">dm</entry>
        <entry key="cm">cm</entry>
        <entry key="millimeters">mm</entry>
        <entry key="micrometers">μm</entry>
        <entry key="nanometers">nm</entry>
        <entry key="picometers">pm</entry>
        <entry key="femtometers">fm</entry>
        <entry key="km">km</entry>
        <entry key="angstrom2">A2</entry>
        <entry key="microm2">μm2</entry>
        <entry key="mm2">mm2</entry>
        <entry key="cm2">cm2</entry>
        <entry key="km2">km2</entry>
        <entry key="ha">ha</entry>
        <entry key="angstrom3">Å³</entry>
        <entry key="L">L</entry>
        <entry key="cm3">cm3</entry>
        <entry key="mm3">mm3</entry>
        <entry key="microm3">μm3</entry>
        <entry key="ppmv">ppmV</entry>
        <entry key="ppbv">ppbV</entry>
        <entry key="pptv">pptV</entry>
        <entry key="cm3.m3">cm3.m3</entry>
        <entry key="microm3.m3">μm3.m3</entry>
        <entry key="megagramsPerCubicMetre">Mg.m-3</entry>
        <entry key="ng.m-3">ng.m-3</entry>
        <entry key="ng.cm-3">ng.cm-3</entry>
        <entry key="g.cm-3">g.cm-3</entry>
        <entry key="g.m-3">g.m-3</entry>
        <entry key="electronsPerCubicAngstrom">eA-3</entry>
        <entry key="DU">DU</entry>
        <entry key="mmmol.m-2">mmmol.m-2</entry>
        <entry key="micromole.mol-1">μmol.mol-1</entry>
        <entry key="nmol.mol-1">nmol.mol-1</entry>
        <entry key="pmol.mol-1">pmol.mol-1</entry>
        <entry key="dn_dr">dn/dr</entry>
        <entry key="pmol">pmol</entry>
        <entry key="nmol">nmol</entry>
        <entry key="micromol">μmol</entry>
        <entry key="millimol">mmol</entry>
        <entry key="julianYear">a_j</entry>
        <entry key="tropicalYear">a_t</entry>
        <entry key="wk">Wk</entry>
        <entry key="day">day</entry>
        <entry key="hr">h</entry>
        <entry key="min">min</entry>
        <entry key="milliseconds">ms</entry>
        <entry key="microseconds">μs</entry>
        <entry key="nanoseconds">ns</entry>
        <entry key="picoseconds">ps</entry>
        <entry key="femtoseconds">fs</entry>
        <entry key="day-1">day-1</entry>
        <entry key="hr-1">h-1</entry>
        <entry key="min-1">min-1</entry>
        <entry key="kHz">kHz</entry>
        <entry key="MHz">MHz</entry>
        <entry key="GHz">GHz</entry>
        <entry key="THz">THz</entry>
        <entry key="molec.cm-3.s-1">molecule.cm-3.s-1</entry>
        <entry key="molec.m-3.s-1">molecule.m-3.s-1</entry>
        <entry key="nmol.cm-3.s-1">nmol.cm-3.s-1</entry>
        <entry key="molec-1.cm3.s-1">molecule-1.cm3.s-1</entry>
        <entry key="molec-1.m3.s-1">molecule-1.m3.s-1</entry>
        <entry key="molec-2.cm6.s-1">molecule-2.cm6.s-1</entry>
        <entry key="molec-2.m6.s-1">molecule-2.m6.s-1</entry>
        <entry key="molec.m-2.s-1">molecule.m-2.s-1</entry>
        <entry key="photon.m-2.s-1">photon.m-2.s-1</entry>
        <entry key="particle.m-2.s-1">particle.m-2.s-1</entry>
        <entry key="mg.m-2.h-1">mg.m-2.h-1</entry>
        <entry key="mg.m-2.d-1">mg.m-2.d-1</entry>
        <entry key="microg.m-2.s-1">μg.m-2.s-1</entry>
        <entry key="cm2.molecule-1">cm2.molecule-1</entry>
        <entry key="m2.molecule-1">m2.molecule-1</entry>
        <entry key="molecule.photon-1">molecule.photon-1</entry>
        <entry key="milliamperes">mA</entry>
        <entry key="kiloVolts">kV</entry>
        <entry key="elementaryCharge">e</entry>
        <entry key="kJ.mol-1.nm-1">kJ.mol-1.nm-1</entry>
        <entry key="electronvolt">eV</entry>
        <entry key="kj.mol-1">kJ.mol-1</entry>
        <entry key="hartree">Eh</entry>
        <entry key="kilowatt">kW</entry>
        <entry key="kPascal">kPa</entry>
        <entry key="hpascal">hPa</entry>
        <entry key="atm">atm</entry>
        <entry key="bar">bar</entry>
        <entry key="degreesC">deg</entry>
        <entry key="reciprocalAngstrom">A-1</entry>
        <entry key="reciprocalMillimeters">mm-1</entry>
        <entry key="cm-1">cm-1</entry>
        <entry key="Petag">Pg</entry>
        <entry key="Tg">Tg</entry>
        <entry key="Mg">Mg</entry>
        <entry key="g">g</entry>
        <entry key="mg">mg</entry>
        <entry key="microgram">μg</entry>
        <entry key="ng">ng</entry>
        <entry key="pg">pg</entry>
        <entry key="dalton">Da</entry>
    </xsl:variable>
    <xsl:variable name="nonsi2UnitsSymbolMap">
        <entry key="kcal.mol-1">kcal.mol-1</entry>
        <entry key="ppm">ppm</entry>
        <entry key="cal.mol-1.K-1">cal.mol-1.K-1</entry>
        <entry key="1e-40.esu2.cm2">1e-40 esu2 cm2</entry>
        <entry key="km.mol-1">km/mole</entry>
        <entry key="au">A.U.</entry>
        <entry key="debye">D</entry>
        <entry key="rydberg">Ry</entry>
        <entry key="bohrmag.cell-1">Bohr magneton / cell</entry>
        <entry key="lattice.constant">alat</entry>
        <entry key="ev.angstrom-1">eV/Å</entry>
        <entry key="meV">meV</entry>
        <entry key="picogram">pm</entry>
        <entry key="amu">amu</entry>
        <entry key="attogram">ag</entry>
        <entry key="angstrom.femtoseconds-1">Å/fs</entry>
        <entry key="angstrom.picoseconds-1">Å/ps</entry>
        <entry key="erg">erg</entry>
        <entry key="picogram-micrometers2.microseconds-2">pg-μm2/μs2</entry>
        <entry key="attogram-nanometers2.nanoseconds-2">ag-nm2/ns2</entry>
        <entry key="bohr.amu-1">bohr/amu</entry>
        <entry key="micrometers.microseconds-1">μm/μs</entry>
        <entry key="nanometers.nanoseconds-1">nm/ns</entry>
        <entry key="picogram.micrometers-microsecond-2">pg/μm-μs2</entry>
        <entry key="picogram.micrometers-microsecond-1">pg/μm-μs</entry>
        <entry key="picocoulombs">pC</entry>
        <entry key="attogram-nanometer.nanosecond-2">ag-nm/ns2</entry>
        <entry key="attogram-nanometer2.nanosecond-2">ag-nm2/ns2</entry>
        <entry key="picogram-micrometers.microseconds-2">pg-μm/μs2</entry>
        <entry key="attogram.nanometers-nanoseconds-2">ag/nm-ns2</entry>
        <entry key="attogram.nanometers-nanoseconds-1">ag/nm-ns</entry>
        <entry key="poise">P</entry>
        <entry key="pascal.s">Pa*s</entry>
        <entry key="statcoulomb">statC</entry>
        <entry key="elementaryCharge.angstrom-1">e/Å</entry>
        <entry key="elementaryCharge.angstrom">e.Å</entry>
        <entry key="picocoulombs-micrometers-1">pC/μm</entry>
        <entry key="elementaryCharge-nanometer">e-nm</entry>
        <entry key="volt.angstrom-1">V/Å</entry>
        <entry key="coulombs.m">C.m</entry>
        <entry key="volt.meter-1">V/m</entry>
        <entry key="volt.cm-1">V/cm</entry>
        <entry key="volt.micrometer-1">V/μm</entry>
        <entry key="volt.nanometer-1">V/nm</entry>
        <entry key="dyne">dyn</entry>
        <entry key="dyne-cm">dyn-cm</entry>
        <entry key="hartree.bohr-1">Eh/Bohr</entry>
        <entry key="electronvolt.atom-1">eV/atom</entry>
        <entry key="electronvolt.angstrom-1">eV/Å</entry>
        <entry key="newton-m">N-m</entry>
        <entry key="kcal.mol-angstrom-1">kcal.mol-angstrom-1</entry>
        <entry key="picogram.micrometers-3">pg/μm3</entry>
        <entry key="attogram.nanometers-3">ag/nm3</entry>
        <entry key="hbar.2-1">hbar/2</entry>
        <entry key="gigapascal">GPa</entry>
        <entry key="electronvolt.angstrom-2">eV/Å2</entry>
        <entry key="debye.angstrom-1.2.amu">(D/A)**2/amu</entry>
        <entry key="step">step(s)</entry>
        <entry key="kilobyte">kB</entry>
        <entry key="cycle">cycle(s)</entry>
        <entry key="angstrom4.amu-1">Å^4/amu</entry>
        <entry key="kJ.mol-1.K-1">kJ/mol/K</entry>
        <entry key="hartree.K-1">Eh/K</entry>
    </xsl:variable>

    <xsl:variable name="atomType2AtomicNumber">
        <entry key="H">1</entry>
        <entry key="He">2</entry>
        <entry key="Li">3</entry>
        <entry key="Be">4</entry>
        <entry key="B">5</entry>
        <entry key="C">6</entry>
        <entry key="N">7</entry>
        <entry key="O">8</entry>
        <entry key="F">9</entry>
        <entry key="Ne">10</entry>
        <entry key="Na">11</entry>
        <entry key="Mg">12</entry>
        <entry key="Al">13</entry>
        <entry key="Si">14</entry>
        <entry key="P">15</entry>
        <entry key="S">16</entry>
        <entry key="Cl">17</entry>
        <entry key="Ar">18</entry>
        <entry key="K">19</entry>
        <entry key="Ca">20</entry>
        <entry key="Sc">21</entry>
        <entry key="Ti">22</entry>
        <entry key="V">23</entry>
        <entry key="Cr">24</entry>
        <entry key="Mn">25</entry>
        <entry key="Fe">26</entry>
        <entry key="Co">27</entry>
        <entry key="Ni">28</entry>
        <entry key="Cu">29</entry>
        <entry key="Zn">30</entry>
        <entry key="Ga">31</entry>
        <entry key="Ge">32</entry>
        <entry key="As">33</entry>
        <entry key="Se">34</entry>
        <entry key="Br">35</entry>
        <entry key="Kr">36</entry>
        <entry key="Rb">37</entry>
        <entry key="Sr">38</entry>
        <entry key="Y">39</entry>
        <entry key="Zr">40</entry>
        <entry key="Nb">41</entry>
        <entry key="Mo">42</entry>
        <entry key="Tc">43</entry>
        <entry key="Ru">44</entry>
        <entry key="Rh">45</entry>
        <entry key="Pd">46</entry>
        <entry key="Ag">47</entry>
        <entry key="Cd">48</entry>
        <entry key="In">49</entry>
        <entry key="Sn">50</entry>
        <entry key="Sb">51</entry>
        <entry key="Te">52</entry>
        <entry key="I">53</entry>
        <entry key="Xe">54</entry>
        <entry key="Cs">55</entry>
        <entry key="Ba">56</entry>
        <entry key="La">57</entry>
        <entry key="Ce">58</entry>
        <entry key="Pr">59</entry>
        <entry key="Nd">60</entry>
        <entry key="Pm">61</entry>
        <entry key="Sm">62</entry>
        <entry key="Eu">63</entry>
        <entry key="Gd">64</entry>
        <entry key="Tb">65</entry>
        <entry key="Dy">66</entry>
        <entry key="Ho">67</entry>
        <entry key="Er">68</entry>
        <entry key="Tm">69</entry>
        <entry key="Yb">70</entry>
        <entry key="Lu">71</entry>
        <entry key="Hf">72</entry>
        <entry key="Ta">73</entry>
        <entry key="W">74</entry>
        <entry key="Re">75</entry>
        <entry key="Os">76</entry>
        <entry key="Ir">77</entry>
        <entry key="Pt">78</entry>
        <entry key="Au">79</entry>
        <entry key="Hg">80</entry>
        <entry key="Tl">81</entry>
        <entry key="Pb">82</entry>
        <entry key="Bi">83</entry>
        <entry key="Po">84</entry>
        <entry key="At">85</entry>
        <entry key="Rn">86</entry>
        <entry key="Fr">87</entry>
        <entry key="Ra">88</entry>
        <entry key="Ac">89</entry>
        <entry key="Th">90</entry>
        <entry key="Pa">91</entry>
        <entry key="U">92</entry>
        <entry key="Np">93</entry>
        <entry key="Pu">94</entry>
        <entry key="Am">95</entry>
        <entry key="Cm">96</entry>
        <entry key="Bk">97</entry>
        <entry key="Cf">98</entry>
        <entry key="Es">99</entry>
        <entry key="Fm">100</entry>
        <entry key="Md">101</entry>
        <entry key="No">102</entry>
        <entry key="Lr">103</entry>
        <entry key="Rf">104</entry>
        <entry key="Db">105</entry>
        <entry key="Sg">106</entry>
        <entry key="Bh">107</entry>
        <entry key="Hs">108</entry>
        <entry key="Mt">109</entry>
        <entry key="Ds">110</entry>
        <entry key="Rg">111</entry>
        <entry key="Cn">112</entry>
        <entry key="Nh">113</entry>
        <entry key="Fl">114</entry>
        <entry key="Mc">115</entry>
        <entry key="Lv">116</entry>
        <entry key="Ts">117</entry>
        <entry key="Og">118</entry>
    </xsl:variable>

    <xsl:variable name="atomTypeOrderedCHFirst">
        <entry key="C"/>
        <entry key="H"/>
        <entry key="Ac"/>
        <entry key="Ag"/>
        <entry key="Al"/>
        <entry key="Am"/>
        <entry key="Ar"/>
        <entry key="As"/>
        <entry key="At"/>
        <entry key="Au"/>
        <entry key="B"/>
        <entry key="Ba"/>
        <entry key="Be"/>
        <entry key="Bh"/>
        <entry key="Bi"/>
        <entry key="Bk"/>
        <entry key="Br"/>
        <entry key="Ca"/>
        <entry key="Cd"/>
        <entry key="Ce"/>
        <entry key="Cf"/>
        <entry key="Cl"/>
        <entry key="Cm"/>
        <entry key="Cn"/>
        <entry key="Co"/>
        <entry key="Cr"/>
        <entry key="Cs"/>
        <entry key="Cu"/>
        <entry key="Db"/>
        <entry key="Ds"/>
        <entry key="Dy"/>
        <entry key="Er"/>
        <entry key="Es"/>
        <entry key="Eu"/>
        <entry key="F"/>
        <entry key="Fe"/>
        <entry key="Fl"/>
        <entry key="Fm"/>
        <entry key="Fr"/>
        <entry key="Ga"/>
        <entry key="Gd"/>
        <entry key="Ge"/>
        <entry key="He"/>
        <entry key="Hf"/>
        <entry key="Hg"/>
        <entry key="Ho"/>
        <entry key="Hs"/>
        <entry key="I"/>
        <entry key="In"/>
        <entry key="Ir"/>
        <entry key="K"/>
        <entry key="Kr"/>
        <entry key="La"/>
        <entry key="Li"/>
        <entry key="Lr"/>
        <entry key="Lu"/>
        <entry key="Lv"/>
        <entry key="Mc"/>
        <entry key="Md"/>
        <entry key="Mg"/>
        <entry key="Mn"/>
        <entry key="Mo"/>
        <entry key="Mt"/>
        <entry key="N"/>
        <entry key="Na"/>
        <entry key="Nb"/>
        <entry key="Nd"/>
        <entry key="Ne"/>
        <entry key="Nh"/>
        <entry key="Ni"/>
        <entry key="No"/>
        <entry key="Np"/>
        <entry key="O"/>
        <entry key="Og"/>
        <entry key="Os"/>
        <entry key="P"/>
        <entry key="Pa"/>
        <entry key="Pb"/>
        <entry key="Pd"/>
        <entry key="Pm"/>
        <entry key="Po"/>
        <entry key="Pr"/>
        <entry key="Pt"/>
        <entry key="Pu"/>
        <entry key="Ra"/>
        <entry key="Rb"/>
        <entry key="Re"/>
        <entry key="Rf"/>
        <entry key="Rg"/>
        <entry key="Rh"/>
        <entry key="Rn"/>
        <entry key="Ru"/>
        <entry key="S"/>
        <entry key="Sb"/>
        <entry key="Sc"/>
        <entry key="Se"/>
        <entry key="Sg"/>
        <entry key="Si"/>
        <entry key="Sm"/>
        <entry key="Sn"/>
        <entry key="Sr"/>
        <entry key="Ta"/>
        <entry key="Tb"/>
        <entry key="Tc"/>
        <entry key="Te"/>
        <entry key="Th"/>
        <entry key="Ti"/>
        <entry key="Tl"/>
        <entry key="Tm"/>
        <entry key="Ts"/>
        <entry key="U"/>
        <entry key="V"/>
        <entry key="W"/>
        <entry key="Xe"/>
        <entry key="Y"/>
        <entry key="Yb"/>
        <entry key="Zn"/>
        <entry key="Zr"/>
    </xsl:variable>
    <!-- Helper functions -->
    <!-- Sqrt function -->
    <xsl:function name="ckbk:sqrt" as="xs:double">
        <xsl:param name="number" as="xs:double"/>
        <xsl:variable name="try"
            select="if ($number lt 100.0) then 1.0
            else if ($number gt 100.0 and $number lt
            1000.0) then 10.0
            else if ($number gt 1000.0 and $number lt
            10000.0) then 31.0
            else 100.00"
            as="xs:decimal"/>

        <xsl:sequence select="if ($number ge 0) then ckbk:sqrt($number,$try,1,20) else -1"/>
    </xsl:function>
    <xsl:function name="ckbk:sqrt" as="xs:double">
        <xsl:param name="number" as="xs:double"/>
        <xsl:param name="try" as="xs:double"/>
        <xsl:param name="iter" as="xs:integer"/>
        <xsl:param name="maxiter" as="xs:integer"/>

        <xsl:variable name="result" select="$try * $try" as="xs:double"/>
        <xsl:sequence
            select="if ($result eq $number or $iter gt $maxiter)
            then $try
            else ckbk:sqrt($number, ($try - (($result - $number)
            div (2 * $try))), $iter + 1, $maxiter)"
        />
    </xsl:function>

    <!-- Simple map example from http://stackoverflow.com/questions/3626118/xslt-creating-a-map-in-xslt -->
    <xsl:function name="helper:printUnitSymbol">
        <xsl:param name="unitQName" as="xs:string?"/>
        <xsl:choose>
            <xsl:when test="not(exists($unitQName)) or compare($unitQName,'') = 0 ">
                <xsl:sequence select="''"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:variable name="namespace" select="substring-before($unitQName,':')"/>
                <xsl:variable name="name" select="substring-after($unitQName,':')"/>
                <xsl:choose>
                    <xsl:when test="matches($namespace,$siUnitsPrefix)">
                        <xsl:sequence select="$siUnitsSymbolMap/entry[@key=$name]"></xsl:sequence>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:choose>
                            <xsl:when test="matches($namespace,$nonsiUnitsPrefix)">
                                <xsl:sequence select="$nonsiUnitsSymbolMap/entry[@key=$name]"></xsl:sequence>
                            </xsl:when>
                            <xsl:otherwise>
                                <xsl:choose>
                                    <xsl:when test="matches($namespace,$nonsi2UnitsPrefix)">
                                        <xsl:sequence select="$nonsi2UnitsSymbolMap/entry[@key=$name]"></xsl:sequence>
                                    </xsl:when>
                                </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:function>

    <!-- Convert atom types to it's correspoding atomic numbers -->
    <xsl:function name="helper:atomType2Number">
        <xsl:param name="atomType" as="xs:string"/>
        <xsl:sequence select="$atomType2AtomicNumber/entry[@key=$atomType]"/>
    </xsl:function>

    <!-- Given a cml molecule, it calculates it's formula on Hill notation-->
    <xsl:function name="helper:calculateHillNotationFormula" >
        <xsl:param name="atomArray"/>
        <xsl:for-each select="$atomTypeOrderedCHFirst/entry">
            <xsl:variable name="element" select="./@key"/>
            <xsl:variable name="elementCount" select="count($atomArray[@elementType=$element])"/>
            <xsl:choose>
                <xsl:when test="$elementCount = 1">
                    <xsl:value-of select="$element"/>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:if test="$elementCount > 1">
                        <xsl:value-of select="concat($element,$elementCount)"/>
                    </xsl:if>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:for-each>
    </xsl:function>

    <xsl:function name="helper:normalizeHillNotation">
        <xsl:param name="formula"/>
        <xsl:variable name="formulaNoWhites" select="replace($formula,' ','')"/>
        <xsl:variable name="formulaNoMiddleOnes"  select="replace($formulaNoWhites, '([A-Za-z]+)1([A-Za-z]+)','$1$2')"/>
        <xsl:value-of select="replace($formulaNoMiddleOnes, '([A-Za-z]+)1$','$1')"/>
    </xsl:function>

    <xsl:function name="helper:isTransitionMetal">
        <xsl:param name="elementNumber"/>
        <xsl:value-of select="( $elementNumber >= 21 and $elementNumber &lt;= 30) or ( $elementNumber >= 39 and $elementNumber &lt;= 48) or ( $elementNumber >= 57 and $elementNumber &lt;= 80)"/>
    </xsl:function>

    <!-- Trim whitespaces from string -->
    <xsl:function name="helper:trim" as="xs:string">
        <xsl:param name="arg" as="xs:string?"/>
        <xsl:sequence select="replace(replace($arg,'\s+$',''),'^\s+','')"/>
    </xsl:function>

    <!-- Given a cml matrix element, it generates a JSON array string -->
    <xsl:function name="helper:matrixToJSON">
        <xsl:param name="matrix"/>
        <xsl:variable name="matrixValues" select="tokenize(helper:trim($matrix/text()),'\s+')"/>
        <xsl:variable name="cols" select="$matrix/@cols | $matrix/@columns"/>
        <xsl:variable name="rows" select="$matrix/@rows"/>
        <xsl:for-each select="1 to $rows">
            <xsl:variable name="outerIndex" select="."/>
            [<xsl:for-each select="1 to $cols">
                <xsl:variable name="innerIndex" select="."/>
                '<xsl:value-of select="$matrixValues[ ($outerIndex - 1) * $cols + ($innerIndex ) ]"/>'<xsl:if test="$innerIndex &lt; $cols">,</xsl:if>
            </xsl:for-each>]<xsl:if test="$outerIndex &lt; $rows">,</xsl:if>
        </xsl:for-each>
    </xsl:function>

    <!-- Given a cml matrix element, it generates a JSON array string, with provided delimiter -->
    <xsl:function name="helper:matrixToJSON">
        <xsl:param name="matrix"/>
        <xsl:param name="delimiter"/>
        <xsl:variable name="matrixValues" select="tokenize(helper:trim($matrix/text()), $delimiter)"/>
        <xsl:variable name="cols" select="$matrix/@cols | $matrix/@columns"/>
        <xsl:variable name="rows" select="$matrix/@rows"/>
        <xsl:for-each select="1 to $rows">
            <xsl:variable name="outerIndex" select="."/>[<xsl:for-each select="1 to $cols"><xsl:variable name="innerIndex" select="."/>'<xsl:value-of select="$matrixValues[ ($outerIndex - 1) * $cols + ($innerIndex ) ]"/>'<xsl:if test="$innerIndex &lt; $cols">,</xsl:if></xsl:for-each>]<xsl:if test="$outerIndex &lt; $rows">,</xsl:if>
        </xsl:for-each>
    </xsl:function>


    <!-- Given a cml matrix element, it generates a JSON array, ordered by columns-->
    <xsl:function name="helper:matrixToJSON2">
        <xsl:param name="matrix"/>
        <xsl:variable name="matrixValues" select="tokenize(helper:trim($matrix/text()),'\s+')"/>
        <xsl:variable name="cols" select="$matrix/@cols | $matrix/@columns"/>
        <xsl:variable name="rows" select="$matrix/@rows"/>
        [<xsl:for-each select="1 to $cols"><xsl:variable name="outerIndex" select="."/>[<xsl:for-each select="$matrixValues[position() mod $cols = $outerIndex mod $cols]"><xsl:variable name="innerIndex" select="position()"/><xsl:value-of select="if(number(.) = 0) then '0' else ."/><xsl:if test="$innerIndex &lt; $rows">,</xsl:if></xsl:for-each>]<xsl:if test="$outerIndex &lt; $cols">,</xsl:if></xsl:for-each>]
    </xsl:function>

    <xsl:function name="helper:correctAtomTypeCase">
        <xsl:param name="atomType"/>
        <xsl:variable name="atomTypeTrimmed" select="helper:trim($atomType)"/>
        <xsl:variable name="correctedAtom">
            <xsl:value-of select="upper-case(substring($atomTypeTrimmed,1,1))"/>
            <xsl:if test="string-length($atomTypeTrimmed) = 2">
                <xsl:value-of select="lower-case(substring($atomTypeTrimmed, 2,1))"/>
            </xsl:if>
            <xsl:if test="string-length($atomTypeTrimmed) = 3">
                <xsl:value-of select="lower-case(substring($atomTypeTrimmed, 2,2))"/>
            </xsl:if>
        </xsl:variable>
        <xsl:sequence select="$correctedAtom"/>
    </xsl:function>

    <xsl:function name="helper:formatScientific">
        <xsl:param name="number"/>
        <xsl:param name="numberOfDecimals"/>
        <!-- 1 to -1 range -->
        <xsl:if test="xs:decimal($number) &lt; 1.0 and xs:decimal($number) &gt; -1.0 ">
            <xsl:variable name="numberStr" select="string(xs:decimal($number))"/>
            <xsl:variable name="isNegative" select="contains($numberStr,'-0.')"/>
            <xsl:variable name="exponent" select="string-length((tokenize(replace($numberStr,'^-?0\.',''),'[1-9]'))[1]) + 1"/>
            <xsl:variable name="base" select="replace($numberStr,'^-?0\.0*','')"/>
            <xsl:variable name="scientificNotation">
                <xsl:choose>
                    <xsl:when test="xs:decimal($number) = 0">
                        0.00
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:if test="$isNegative">-</xsl:if><xsl:value-of select="concat(substring($base,1,1), '.', substring($base,2,$numberOfDecimals))"/>E-<xsl:value-of select="format-number(number($exponent),'00')"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:variable>
            <xsl:value-of select="$scientificNotation"/>
        </xsl:if>
        <!-- greater than 1 to less than -1 range -->
        <xsl:if test="xs:decimal($number) &gt;= 1.0 or xs:decimal($number) &lt;= -1.0 ">
            <xsl:variable name="numberStr" select="string(xs:decimal($number))"/>
            <xsl:variable name="isNegative" select="matches($numberStr,'\-[0-9]+.[0-9]+')"/>
            <xsl:variable name="exponent" select="string-length(replace(replace($numberStr,'\.[0-9]+',''),'-','')) -1"/> <!-- Remove decimal part and count integer positions -->

            <xsl:variable name="base" select="replace($numberStr,'[\.-]','')"/>
            <xsl:variable name="scientificNotation">
                <xsl:if test="$isNegative">-</xsl:if>
                <xsl:value-of select="concat(substring($base,1,1), '.', substring($base,2,$numberOfDecimals))"/><xsl:if test="$exponent!= 0">E<xsl:value-of select="format-number(number($exponent),'00')"/></xsl:if>
            </xsl:variable>
            <xsl:value-of select="$scientificNotation"/>
        </xsl:if>
    </xsl:function>

    <!-- Helper mathematical functions  -->
    <!-- <xsl:include href="../math/fxsl/sqrt.xsl"/> -->
    <!-- <xsl:include href="../math/fxsl/arcTrignm.xsl"/> -->
    <!-- <xsl:include href="../math/fxsl/transform-and-sum.xsl"/> -->

    <!-- Given a vector, this function calculates it's length -->
    <xsl:function name="helper:calcAxisLength">
        <xsl:param name="a"/>
        <xsl:value-of select="format-number(ckbk:sqrt(number($a[1])*number($a[1]) + number($a[2])*number($a[2]) + number($a[3])*number($a[3])),'#0.00000000')"/>
    </xsl:function>

    <!-- Given two vectors, this function returns the angle between them in degrees (rounded to two decimals) -->
<!--    <xsl:function name="helper:calcAxesAngle">
        <xsl:param name="a"/>
        <xsl:param name="b"/>
        <xsl:variable name="aAxisLength" select="helper:calcAxisLength($a)"/>
        <xsl:variable name="bAxisLength" select="helper:calcAxisLength($b)"/>

        <xsl:variable name="result">
            <xsl:call-template name="arccos">
                <xsl:with-param name="pX" select="((number($a[1]) * number($b[1])) + (number($a[2]) * number($b[2])) + (number($a[3]) * number($b[3]))) div ($aAxisLength * $bAxisLength)"/>
            </xsl:call-template>
        </xsl:variable>
        <xsl:variable name="resultGrad" select="57.2957795 * $result"/>
        <xsl:value-of select="format-number( round($resultGrad*100) div 100 ,'#0.00')"/>

    </xsl:function>-->

    <!-- Convert fractionals coordinates to cartesian, using unit cell vectors and atom fractional atom coordinate -->
    <xsl:function name="helper:fractionalToCartesian" xml:space="default">
        <xsl:param name="a"/>
        <xsl:param name="b"/>
        <xsl:param name="c"/>
        <xsl:param name="frac"/>
        <xsl:element name="coords" xml:space="default">
            <xsl:attribute name="x3" select="format-number((number($frac[1]) * number($a[1])) + (number($frac[2]) * number($b[1])) + (number($frac[3]) * number($c[1])),'#0.00000000')"/>
            <xsl:attribute name="y3" select="format-number((number($frac[1]) * number($a[2])) + (number($frac[2]) * number($b[2])) + (number($frac[3]) * number($c[2])),'#0.00000000')"/>
            <xsl:attribute name="z3" select="format-number((number($frac[1]) * number($a[3])) + (number($frac[2]) * number($b[3])) + (number($frac[3]) * number($c[3])),'#0.00000000')"/>
        </xsl:element>
    </xsl:function>

    <!-- This function reads a group of nodes which contain one or multiple text values, separated by blanks and builds an array element -->
    <xsl:function name="helper:arrayFromStringNodes">
        <xsl:param name="dictRef"/>
        <xsl:param name="dataType"/>
        <xsl:param name="nodes"/>

        <xsl:variable name="values">
            <xsl:for-each select="$nodes">
                <xsl:value-of select="tokenize(helper:trim(.),'\s+')"/><xsl:text> </xsl:text>
            </xsl:for-each>
        </xsl:variable>

        <xsl:element name="array">
            <xsl:attribute name="dataType"><xsl:value-of select="$dataType"/></xsl:attribute>
            <xsl:attribute name="dictRef"><xsl:value-of select="$dictRef"/></xsl:attribute>
            <xsl:attribute name="size"><xsl:value-of select="count(tokenize(helper:trim($values),'\s+'))"/></xsl:attribute>
            <xsl:value-of select="helper:trim($values)"/>
        </xsl:element>
    </xsl:function>

    <xsl:function name="helper:arrayFromStringNodesWithDelimiter">
        <xsl:param name="dictRef"/>
        <xsl:param name="dataType"/>
        <xsl:param name="nodes"/>
        <xsl:param name="delimiter"/>
        <xsl:param name="delimiterRegExp"/>


        <xsl:variable name="values">
            <xsl:for-each select="$nodes">
                <xsl:value-of select="helper:trim(.)"/><xsl:if test="position() != last()"><xsl:value-of select="$delimiter"/></xsl:if>
            </xsl:for-each>
        </xsl:variable>

        <xsl:element name="array">
            <xsl:attribute name="delimiter"><xsl:value-of select="$delimiter"/></xsl:attribute>
            <xsl:attribute name="dataType"><xsl:value-of select="$dataType"/></xsl:attribute>
            <xsl:attribute name="dictRef"><xsl:value-of select="$dictRef"/></xsl:attribute>
            <xsl:attribute name="size"><xsl:value-of select="count(tokenize($values,$delimiterRegExp))"/></xsl:attribute>
            <xsl:value-of select="$values"/>
        </xsl:element>
    </xsl:function>

    <!-- This function reads a group of nodes which contain one or multiple text values, separated by blanks and builds a matrix element -->
    <!-- It also calculates rows and cols, one such fields must be passed to the function, the other will be set to -1 -->
    <xsl:function name="helper:matrixFromStringNodes">
        <xsl:param name="dictRef"/>
        <xsl:param name="dataType"/>
        <xsl:param name="nodes"/>
        <xsl:param name="refRows"/>
        <xsl:param name="refCols"/>

        <xsl:variable name="values">
            <xsl:for-each select="$nodes">
                <xsl:value-of select="tokenize(helper:trim(.),'\s+')"/><xsl:text> </xsl:text>
            </xsl:for-each>
        </xsl:variable>

        <xsl:variable name="rows">
            <xsl:choose>
                <xsl:when test="$refRows != -1">
                    <xsl:value-of select="$refRows"/>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:value-of select="count(tokenize(helper:trim($values),'\s+')) div $refCols"/>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:variable>

        <xsl:variable name="cols">
            <xsl:choose>
                <xsl:when test="$refCols != -1">
                    <xsl:value-of select="$refCols"/>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:value-of select="count(tokenize(helper:trim($values),'\s+')) div $refRows"/>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        <xsl:element name="matrix">
            <xsl:attribute name="dataType"><xsl:value-of select="$dataType"/></xsl:attribute>
            <xsl:attribute name="dictRef"><xsl:value-of select="$dictRef"/></xsl:attribute>
            <xsl:attribute name="rows"><xsl:value-of select="$rows"/></xsl:attribute>
            <xsl:attribute name="cols"><xsl:value-of select="$cols"/></xsl:attribute>
            <xsl:value-of select="helper:trim($values)"/>
        </xsl:element>

    </xsl:function>
    <xsl:function name="helper:readBoolean">
        <xsl:param name="value"/>
        <xsl:choose>
            <xsl:when test="compare(helper:trim($value), 'T') = 0">
                <xsl:sequence select="true()"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:sequence select="false()"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:function>

    <xsl:function name="helper:countNumberOfValues">
          <xsl:param name="value"/>
        <xsl:value-of select="count(tokenize(helper:trim($value),'\s+'))"/>
    </xsl:function>

    <xsl:function name="helper:sumStringArray">
        <xsl:param name="values"/>
        <xsl:variable name="sum">
            <xsl:for-each select="1 to count($values)">
                <xsl:variable name="outerIndex" select="."/>
                <xsl:element name="value">
                    <xsl:value-of select="$values[$outerIndex]"/>
                </xsl:element>
            </xsl:for-each>
        </xsl:variable>
        <xsl:value-of select="sum($sum/*:value)"/>
    </xsl:function>

    <!-- Conversion tables taken from http://www.colby.edu/chemistry/PChem/Hartree.html -->
    <xsl:variable name="unitsKcalMol" select="'nonsi2:kcal.mol-1'"/>
    <xsl:variable name="unitsKJMol"   select="'nonsi:kj.mol-1'"/>
    <xsl:variable name="unitsEV"      select="'nonsi:electronvolt'"/>
    <xsl:variable name="unitsHartree" select="'nonsi:hartree'"/>
    <xsl:variable name="unitsCM1"     select="'nonsi:cm-1'"/>

    <xsl:variable name="auToKcalMol" select="627.51"/>
    <xsl:variable name="auToKJMol"   select="2625.50"/>
    <xsl:variable name="auToEV"      select="27.211"/>
    <xsl:variable name="auToCm1"     select="219470"/>

    <xsl:variable name="evToKcalMol" select="23.061"/>
    <xsl:variable name="evToKJMol"   select="96.485"/>
    <xsl:variable name="evToHartree" select="0.036749"/>
    <xsl:variable name="evToCm1"     select="8065.5"/>

    <xsl:variable name="cm1ToKcalMol" select="0.0028591"/>
    <xsl:variable name="cm1ToKJMol"   select="0.011963"/>
    <xsl:variable name="cm1ToHartree" select="0.0000045563"/>
    <xsl:variable name="cm1ToEV"      select="0.00012398"/>

    <xsl:variable name="kcalMolToCm1" select="349.76"/>
    <xsl:variable name="kcalMolToKJMol" select="4.184"/>
    <xsl:variable name="kcalMolToHartree" select="0.0015936"/>
    <xsl:variable name="kcalMolToEV" select="0.043364"/>

    <xsl:variable name="kJMolToCm1" select="83.593"/>
    <xsl:variable name="kJMolToKcalMol" select="0.23901"/>
    <xsl:variable name="kJMolToHartree" select="0.00038088"/>
    <xsl:variable name="kJMolToEV" select="0.010364"/>

    <xsl:function name="helper:convertEnergyUnits">
        <!-- Conversion table taken from http://www.ucl.ac.uk/~uccaati/Energy.html -->
        <xsl:param name="valueNode"/>
        <xsl:param name="toUnits"/>
        <xsl:if test="exists($valueNode)">
            <xsl:choose>
                <xsl:when test="compare($valueNode/@units, $unitsHartree) = 0">
                    <xsl:choose>
                        <xsl:when test="compare($toUnits,$unitsKcalMol) = 0">
                            <xsl:value-of select="number($valueNode) * $auToKcalMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits,$unitsKJMol) = 0">
                            <xsl:value-of select="number($valueNode) * $auToKJMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits,$unitsEV) = 0">
                            <xsl:value-of select="number($valueNode) * $auToEV"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits,$unitsCM1) = 0">
                            <xsl:value-of select="number($valueNode) * $auToCm1"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits,$unitsHartree) = 0">
                            <xsl:value-of select="$valueNode"/>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>

                <xsl:when test="compare($valueNode/@units,$unitsEV) = 0">
                    <xsl:choose>
                        <xsl:when test="compare($toUnits, $unitsKcalMol) = 0">
                            <xsl:value-of select="number($valueNode) * $evToKcalMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsKJMol) = 0">
                            <xsl:value-of select="number($valueNode) * $evToKJMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsHartree) = 0">
                            <xsl:value-of select="number($valueNode) * $evToHartree"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsCM1) = 0">
                            <xsl:value-of select="number($valueNode) * $evToCm1"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsEV) = 0">
                            <xsl:value-of select="$valueNode"/>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>

                <xsl:when test="compare($valueNode/@units,$unitsCM1) = 0">
                    <xsl:choose>
                        <xsl:when test="compare($toUnits, $unitsKcalMol) = 0">
                            <xsl:value-of select="number($valueNode) * $cm1ToKcalMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsKJMol) = 0">
                            <xsl:value-of select="number($valueNode) * $cm1ToKJMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsHartree) = 0">
                            <xsl:value-of select="number($valueNode) * $cm1ToHartree"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsEV) = 0">
                            <xsl:value-of select="number($valueNode) * $cm1ToEV"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsCM1) = 0">
                            <xsl:value-of select="$valueNode"/>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>

                <xsl:when test="compare($valueNode/@units,$unitsKcalMol) = 0">
                    <xsl:choose>
                        <xsl:when test="compare($toUnits, $unitsCM1) = 0">
                            <xsl:value-of select="number($valueNode) * $kcalMolToCm1"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsKJMol) = 0">
                            <xsl:value-of select="number($valueNode) * $kcalMolToKJMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsHartree) = 0">
                            <xsl:value-of select="number($valueNode) * $kcalMolToHartree"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsEV) = 0">
                            <xsl:value-of select="number($valueNode) * $kcalMolToEV"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits,$unitsKcalMol ) = 0">
                            <xsl:value-of select="$valueNode"/>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>

                <xsl:when test="compare($valueNode/@units,$unitsKJMol) = 0">
                    <xsl:choose>
                        <xsl:when test="compare($toUnits, $unitsCM1) = 0">
                            <xsl:value-of select="number($valueNode) * $kJMolToCm1"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsKcalMol ) = 0">
                            <xsl:value-of select="number($valueNode) * $kJMolToKcalMol"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsHartree) = 0">
                            <xsl:value-of select="number($valueNode) * $kJMolToHartree"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsEV) = 0">
                            <xsl:value-of select="number($valueNode) * $kJMolToEV"/>
                        </xsl:when>
                        <xsl:when test="compare($toUnits, $unitsKJMol ) = 0">
                            <xsl:value-of select="$valueNode"/>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>

                <xsl:otherwise>
                    <xsl:value-of><xsl:text>Unhandled energy unit type</xsl:text></xsl:value-of>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:if>
    </xsl:function>



    <xsl:function name="helper:normalizeMethods">
        <xsl:param name="methods" as="xs:string*"/>

        <xsl:variable name="normalizedMethods">
            <xsl:for-each select="$methods">
                <xsl:variable name="method" select="."/>
                <xsl:choose>
                    <xsl:when test="matches(upper-case($method), '^.*PBEC\s*PBEX.*')">
                        <xsl:value-of select="'PBE '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*BECKE88\s*PERDEW86.*')">
                        <xsl:value-of select="'BP86 '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*BECKE88\s*PERDEW86.*')">
                        <xsl:value-of select="'BP86 '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*BECKE88\s*LYP.*') or matches(upper-case($method),'.*BLYP.*')">
                        <xsl:value-of select="'BLYP '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*B3-?LYP\*.*')">
                        <xsl:value-of select="'B3LYP* '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*B3-?LYP.*')">
                        <xsl:value-of select="'B3LYP '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*OPTX\s*PBEC.*')">
                        <xsl:value-of select="'OPBE '"/>
                    </xsl:when>
                    <xsl:when test="matches(upper-case($method),'^.*B-?P.*')">
                        <xsl:value-of select="'BP86 '"/>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:value-of select="concat($method, ' ')"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:for-each>
        </xsl:variable>
        <xsl:value-of select="helper:trim($normalizedMethods)"/>
    </xsl:function>

    <!-- Utility to append an html template with a molden viewer -->
    <xsl:function name="helper:appendMoldenViewerCode">
        <xsl:param name="moldenurl"/>
        <xsl:param name="webrootpath"/>

        <script type="text/javascript">
            $(document).ready(function(){
                $( "#viewMoldenDiv" ).draggable();
            });
            var moldenOrbitalUrl = '<xsl:value-of select="$moldenurl"/>';
            var j2sPath = "<xsl:value-of select="$webrootpath"/>/xslt/jsmol/j2s";
        </script>
        <div id="viewMoldenDiv">
            <xsl:copy-of select="document('../../cml2htmlViewMolden.html')" />
        </div>
    </xsl:function>


</xsl:stylesheet>
