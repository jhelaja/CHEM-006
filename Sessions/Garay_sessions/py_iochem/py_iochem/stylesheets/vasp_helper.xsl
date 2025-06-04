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
    xmlns:vasp="https://www.vasp.at/"
    exclude-result-prefixes="xs"
    version="2.0">
    
    <!-- Calculation type related variabes -->
    <xsl:variable name="vasp:GeometryOptimization" select="'Geometry optimization'" />
    <xsl:variable name="vasp:SinglePoint" select="'Single point'" />
    <xsl:variable name="vasp:MolecularDynamics" select="'Ab-Initio Molecular Dynamics'" />
    <xsl:variable name="vasp:FrequencyCalculus" select="'Frequencies'" />
    <xsl:variable name="vasp:ImprovedDimerMethod" select="'Improved Dimer Method'" />
    <xsl:variable name="vasp:NudgedElasticBand" select="'Nudged Elastic Band'"/>
    <xsl:variable name="vasp:NotAvailable" select="'N/A'" />
    
    <xsl:variable name="vasp:ClosedShell" select="'Closed shell'"/>
    <xsl:variable name="vasp:OpenShell" select="'Open shell'"/>
    
    <xsl:function name="vasp:getCalcType">
        <xsl:param name="ibrion"/>
        <xsl:choose>
            <xsl:when test="count($ibrion) > 1 and exists($ibrion[text() = '44'])"><xsl:value-of select="$vasp:ImprovedDimerMethod"/></xsl:when>
            <xsl:when test="count($ibrion) > 1"><xsl:value-of select="$vasp:NudgedElasticBand"/></xsl:when>
            <xsl:otherwise>
                <xsl:choose>
                    <xsl:when test="$ibrion = -1"><xsl:value-of select="$vasp:SinglePoint"/></xsl:when>
                    <xsl:when test="$ibrion = 0"><xsl:value-of select="$vasp:MolecularDynamics"/></xsl:when>
                    <xsl:when test="$ibrion &gt; 0 and $ibrion &lt; 4"><xsl:value-of select="$vasp:GeometryOptimization"/></xsl:when>
                    <xsl:when test="$ibrion &gt; 4 and $ibrion &lt; 9"><xsl:value-of select="$vasp:FrequencyCalculus"/></xsl:when>
                    <xsl:when test="$ibrion = 44"><xsl:value-of select="$vasp:ImprovedDimerMethod"/></xsl:when>
                    <xsl:otherwise><xsl:value-of select="$vasp:NotAvailable"/></xsl:otherwise>
                </xsl:choose>                                
            </xsl:otherwise>
        </xsl:choose>
    </xsl:function>    
    
    <xsl:function name="vasp:getShellType">
        <xsl:param name="ispin"/>
        <xsl:choose>
            <xsl:when test="$ispin = 1"><xsl:value-of select="$vasp:ClosedShell"/></xsl:when>
            <xsl:when test="$ispin = 2"><xsl:value-of select="$vasp:OpenShell"/></xsl:when>            
        </xsl:choose>
    </xsl:function>
    
    <xsl:function name="vasp:getFunctional">
        <xsl:param name="gga"/>
        <xsl:param name="lhfcalc"/>
        
        <xsl:param name="aggac"/>         
        <xsl:param name="hfscreen"/>        
        <xsl:param name="luseVdw"/>
        <xsl:param name="zabVdw"/>
        <xsl:param name="param1"/>
        <xsl:param name="param2"/>        
        <xsl:param name="ldau"/>
        <xsl:param name="aexx"/>
        <xsl:param name="aggax"/>
        <xsl:param name="aldac"/>
        
        <xsl:choose>
            <xsl:when test="compare($lhfcalc,'true')=0">
                <xsl:choose>
                    <xsl:when test="$hfscreen=0.2">HSE06</xsl:when>
                    <xsl:when test="$hfscreen=0.3">HSE03</xsl:when>
                    <xsl:when test="compare($gga,'B3')=0 and $aexx = 0.2 and $aggax = 0.72 and $aggac = 0.81 and $aldac = 0.19">B3LYP</xsl:when>
                    <xsl:when test="$aexx = 1.0 and $aldac = 0.0 and $aggac = 0.0">HF</xsl:when>     
                    <xsl:otherwise>PBE0</xsl:otherwise>
                </xsl:choose> 
                <xsl:if test="exists($aexx)">
                    <xsl:value-of select="concat(' AEXX=', round($aexx * 100), '%')"/>    
                </xsl:if>                    
            </xsl:when>
            <xsl:when test="compare($luseVdw , 'true')=0 and $aggac = 0.0">
                <xsl:choose>
                    <xsl:when test="compare($gga,'RE')=0">vdW-DF</xsl:when>
                    <xsl:when test="compare($gga,'OR')=0">optPBE-vdW</xsl:when>
                    <xsl:when test="compare($gga,'BO')=0 and round($param1 * 1000) div 1000 = 0.183 and round($param2 * 100) div 100 = 0.22">optB88-vdW</xsl:when>
                    <xsl:when test="compare($gga,'MK')=0 and round($param1 * 10000) div 10000 = 0.1234 and $param2 = 1.0">optB86d-vdW</xsl:when>
                    <xsl:when test="compare($gga,'MK')=0 and round($param1 * 10000) div 10000 = 0.1234 and round($param2 * 1000000) div 1000000 = 0.711357">rev-vdW-DF2</xsl:when>
                    <xsl:when test="compare($gga,'ML')=0 and $zabVdw = -1.8867">vdW-DF2</xsl:when>
                    <xsl:otherwise>N/A</xsl:otherwise>
                </xsl:choose>
            </xsl:when> 
            <xsl:when test="compare($gga,'91')=0">PW91</xsl:when>
            <xsl:when test="compare($gga,'PE')=0">PBE</xsl:when>
            <xsl:when test="compare($gga,'RP')=0">rPBE</xsl:when>
            <xsl:when test="compare($gga,'AM')=0">AM05</xsl:when>
            <xsl:when test="compare($gga,'PS')=0">PBEsol</xsl:when>                          
            <xsl:otherwise>N/A</xsl:otherwise>             
        </xsl:choose>
        <xsl:if test="compare($ldau,'true')=0">
            <xsl:text>+U</xsl:text>
        </xsl:if>         
    </xsl:function>   
    
    
</xsl:stylesheet>