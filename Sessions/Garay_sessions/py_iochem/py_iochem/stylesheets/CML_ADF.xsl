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
<xsl:stylesheet  
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
    xmlns:xd="http://www.oxygenxml.com/ns/doc/xsl"
    xmlns:cml="http://www.xml-cml.org/schema"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:ckbk="http://my.safaribooksonline.com/book/xml/0596009747/numbers-and-math/77"
    xmlns:adf="http://www.scm.com/ADF/"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"
    
    xpath-default-namespace="http://www.xml-cml.org/schema" exclude-result-prefixes="xs xd cml ckbk adf helper cmlx"
    version="2.0">
    <xd:doc scope="stylesheet">
        <xd:desc>            
            <xd:p><xd:b>Created on:</xd:b> June 18, 2024</xd:p>
            <xd:p><xd:b>Author:</xd:b>Moisés Álvarez Moreno and Diego Garay Ruiz</xd:p>
            <xd:p><xd:b>Center:</xd:b>Institute of Chemical Research of Catalonia</xd:p>
        </xd:desc>       
    </xd:doc>
    <xsl:include href="adf_helper.xsl"/>
    <xsl:include href="chem_helper.xsl"/>
    <xsl:output method="text" omit-xml-declaration="yes" indent="no" />
    <xsl:strip-space elements="*"/>
    
    <xsl:param name="title"/>
    <xsl:param name="author"/>    
    
    <xsl:variable name="runType" select="concat((//parameter[@dictRef='cc:runtype']/scalar/text())[last()], '')"/>
    <xsl:variable name="hasVibrations" select="exists(//property[@dictRef='cc:frequencies'])" /> 
    <xsl:variable name="isQuild" select="exists(//module[@cmlx:templateRef='quild.iteration'])"/>
    <xsl:variable name="isNMR" select="exists(//module[@cmlx:templateRef='nucleus'])"/>
    <xsl:variable name="calcType" select="adf:getCalcType($runType,$isQuild,$isNMR)"/>
    
    <!-- Environment module -->
    <xsl:variable name="programParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:program']"/>
    <xsl:variable name="versionParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:programVersion']"/>
    <!-- Initializacion module -->
    <xsl:variable name="solvation" select="(//module[@dictRef='cc:initialization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='solvation'])[1]"/>
    <xsl:variable name="fragmentFiles" select="(//module[@dictRef='cc:initialization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='fragment.files'])[1]"/>
    <xsl:variable name="initialMolecule" select="(//module[@dictRef='cc:initialization' and child::molecule])[last()]//molecule"/>
    <xsl:variable name="functional" select="distinct-values(//module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='cc:functional']/scalar)"/>
    <xsl:variable name="parameters" select="//module[@dictRef='cc:initialization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='parameters']"/>
    <xsl:variable name="symmetry" select="distinct-values(//module[@dictRef='cc:calculation']//module[@cmlx:templateRef='symmetry']/scalar[@dictRef='a:symmetry'])"/>   

    <!-- Geometry -->
    <xsl:variable name="finalMolecule" select="((//module[@dictRef='cc:finalization' and child::molecule])[last()]//molecule)[1]"/>
    <!-- Thermochemistry -->
    <xsl:variable name="thermochemistryProperty" select="//module[@id='finalization']/propertyList/property[@dictRef='cc:thermochemistry']"/>
    <xsl:variable name="temperature" select="($thermochemistryProperty[1]//scalar[@dictRef='cc:temp'])[1]"/>
    <xsl:variable name="pressure" select="($thermochemistryProperty[1]//scalar[@dictRef='cc:press'])[1]"/>
    <!-- Finalization module -->
    <xsl:variable name="converged" select="if(exists(//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='logfile']/scalar[@dictRef='a:converged'])[last()]) then 
                                                (//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='logfile']/scalar[@dictRef='a:converged'])[last()]
                                            else
                                                //list[@cmlx:templateRef='step']/scalar[@dictRef='x:converged' and upper-case(text())= 'CONVERGED'] "/>    
    <xsl:variable name="charge" select="if(exists((//module[@id='finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='logfile']/scalar[@dictRef='a:charge'])[last()])) then
                                            (//module[@id='finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='logfile']/scalar[@dictRef='a:charge'])[last()] 
                                        else 
                                            //module[@cmlx:templateRef='symmetry']/scalar[@dictRef='a:charge']"/>
    <xsl:variable name="spinpolarization" select="(//scalar[@dictRef='a:spinPolarization'])[last()]"/>
    
    
    <xsl:variable name="multiplicity" select="if(exists($spinpolarization)) then
                                                    string(number($spinpolarization)+1)
                                              else if(exists(//module[@cmlx:templateRef='logfile'])) then
                                                    '1'
                                              else '' "/>

    <xsl:template match="/">              
        <xsl:call-template name="generalInfo"/>                                                                      
        <xsl:variable name="molecule" select="
            if(exists($finalMolecule)) 
                then $finalMolecule
            else
                $initialMolecule
            "/>                                                             
            <xsl:call-template name="atomicCoordinates">
                <xsl:with-param name="molecule"     select="$molecule"/>
                <xsl:with-param name="fragmentFiles" select="$fragmentFiles"/>                                                                                        
            </xsl:call-template>                                        
            <xsl:variable name="inchi" select=".//cml:module[@dictRef='cc:finalization']/cml:molecule/cml:formula[@convention='iupac:inchi']/@inline"/>
            <xsl:variable name="molecularFormula" select=".//cml:module[@dictRef='cc:finalization']/cml:molecule/cml:formula/@concise"/>
        
            <xsl:value-of select="concat('inchi:',$inchi,'#;#','molecularFormula:',$molecularFormula,'#;#')"/>
            <xsl:call-template name="chargemultiplicity"/>                              
            <xsl:call-template name="solvatationInfo"/>
                   
            <xsl:for-each select="//cml:module[@dictRef='cc:job']">                                
               <xsl:variable name="scfConverged" select="exists(.//module[@id='calculation']//module[@cmlx:templateRef='scf']/scalar[@dictRef='cc:scfConverged'])"/>  
               <xsl:if test="exists($scfConverged)">
                   <small> SCF Converged</small>
               </xsl:if>                                               
               <xsl:call-template name="bondingEnergy"/>
               <xsl:call-template name="fitTest"/>
<!--               <xsl:call-template name="molecularOrbitals"/>-->
               <xsl:call-template name="orbitalEnergy"/>
               <xsl:call-template name="mullikenCharges"/>
               <xsl:call-template name="mdc"/>
               <xsl:call-template name="multipoleMoments"/>
               <xsl:call-template name="s2"/>
               <xsl:call-template name="frequencyIntensities"/>
               <xsl:call-template name="frequencyDisplacements"/>
               <xsl:call-template name="frequencies"/>
               <xsl:call-template name="zeropointenergy"/> 
               <xsl:call-template name="thermochemistry"/>
<!--               <xsl:call-template name="excitations"/>-->
<!--               <xsl:call-template name="rotatorystrengths"/>-->
<!--               <xsl:call-template name="nmr"/>                                   -->
<!--               <xsl:call-template name="timing"/>                                                                                                        -->
            </xsl:for-each>
     <!--       <xsl:call-template name="printLicense"/>  -->                          
    </xsl:template>

    <!-- General Info -->
    <xsl:template name="generalInfo">
            <xsl:if test="$title">
                <xsl:value-of select="concat('title:',$title,'#;#')"/>
            </xsl:if>
            <xsl:value-of select="concat('program:',$programParameter/scalar/text(),'#;#')"/>                                        
            <xsl:value-of select="concat('programVersion:',$versionParameter/scalar/text(),'#;#')" />
            <xsl:if test="$author">
                <xsl:value-of select="concat('author:',$author,'#;#')"/>
            </xsl:if>
            <xsl:value-of select="concat('formula:',(//formula/@concise)[1],'#;#')"/>

            <xsl:value-of select="concat('calcType:',$calcType)"/>
            <xsl:choose>
                <xsl:when test="exists($solvation)"><xsl:text>_solv,#;#</xsl:text></xsl:when>
                <xsl:otherwise><xsl:text>_gas,'#;#'</xsl:text></xsl:otherwise>
            </xsl:choose>                                       
            <xsl:if test="exists($functional)">
                <xsl-text>methods:</xsl-text>
                <xsl:for-each select="$functional">
                    <xsl:value-of select="concat(.,' ')"/>
                </xsl:for-each>
                <xsl:text>#;#</xsl:text>                    
            </xsl:if>                
            <xsl:if test="exists($parameters//scalar[@dictRef='cc:functional'])">
                <xsl-text>methods:</xsl-text>
                <xsl:for-each select="distinct-values($parameters//scalar[@dictRef='cc:functional'])">
                        <xsl:value-of select="."/><xsl:text> </xsl:text>    
                 </xsl:for-each>
                <xsl-text>#;#</xsl-text>
            </xsl:if>
            <xsl:if test="exists($parameters//scalar[@dictRef='a:relcor'])">
                <xsl-text>relativisticCorrections:</xsl-text>
                <xsl:for-each select="distinct-values($parameters//scalar[@dictRef='a:relcor'])">
                    <xsl:value-of select="."/><xsl:text> </xsl:text>    
                </xsl:for-each>
                <xsl-text>#;#</xsl-text>
            </xsl:if>
            <xsl:if test="exists($parameters//scalar[@dictRef='a:coretreat'])">
                <xsl-text>coreTreatment:</xsl-text>
                <xsl:for-each select="distinct-values($parameters//scalar[@dictRef='a:coretreat'])">
                    <xsl:value-of select="."/><xsl:text> </xsl:text>
                </xsl:for-each>     
                <xsl-text>#;#</xsl-text>
            </xsl:if>
            <xsl:if test="exists($parameters//scalar[@dictRef='a:electricField'])">
               <xsl-text>electricField:</xsl-text>
                <xsl:for-each select="distinct-values($parameters//scalar[@dictRef='a:electricField'])">
                    <xsl:value-of select="."/><xsl:text> </xsl:text>    
                </xsl:for-each>
                <xsl-text>#;#</xsl-text>
            </xsl:if>
            <xsl:if test="exists($parameters//scalar[@dictRef='a:zeeman'])">
                <xsl-text>hyperfineInteraction:</xsl-text>
                <xsl:for-each select="distinct-values($parameters//scalar[@dictRef='a:zeeman'])">
                    <xsl:value-of select="."/><xsl:text> </xsl:text>    
                </xsl:for-each>
                <xsl-text>#;#</xsl-text>
            </xsl:if>
            <xsl:if test="exists($symmetry)">
                <xsl-text>symmetry:</xsl-text>
                <xsl:for-each select="$symmetry">
                    <xsl:value-of select="."/><xsl:text> </xsl:text>
                </xsl:for-each>
                <xsl-text>#;#</xsl-text>
            </xsl:if>            
            <xsl:if test="exists($thermochemistryProperty)">                                                                          
                <xsl:value-of select="concat('temperature:',$temperature/text(),'#;#')"></xsl:value-of>   
                <xsl:value-of select="concat('temperatureUnits:',helper:printUnitSymbol($temperature/@units),'#;#')"></xsl:value-of>
                <xsl:value-of select="concat('pressure:',$pressure/text(),'#;#')"/>
                <xsl:value-of select="concat('pressureUnits:',helper:printUnitSymbol($pressure/@units),'#;#')"></xsl:value-of>
            </xsl:if>
    </xsl:template>

    <!-- Atomic coordinates -->
    <xsl:template name="atomicCoordinates">        
        <xsl:param name="molecule"/>
        <xsl:param name="fragmentFiles"/>   
        <xsl:variable name="collapseAccordion" select="if(count($molecule/cml:atomArray/cml:atom) > 10) then '' else 'in'"/>
        
        <xsl:if test="contains($calcType,$adf:GeometryOptimization)">
            <xsl:choose>
                <xsl:when test="matches($converged/text(),$adf:GeometryConverged)">
                    <xsl-text>optConvergeStatus:True#;#</xsl-text>    
                </xsl:when>
                <xsl:otherwise>
                    <xsl:choose>
                        <xsl:when test="matches($converged/text(),$adf:GeometryNotConverged)">                                                
                            <xsl-text>optConvergeStatus:False#;#</xsl-text>      
                        </xsl:when>
                    </xsl:choose>                                        
                </xsl:otherwise>
            </xsl:choose>                                                           
        </xsl:if>                    
                         
        <xsl-text>geometryCartesian:</xsl-text>  
        <xsl:for-each select="$molecule/cml:atomArray/cml:atom">                                                        
            <xsl:variable name="outerIndex" select="position()"/>
            <xsl:variable name="elementType" select="@elementType"/>                                                                                                       
            <xsl:variable name="id" select="@id"/>
            <xsl:value-of select="concat($elementType,' ',format-number(@x3, '#0.0000'),' ',format-number(@y3, '#0.0000'),' ',format-number(@z3, '#0.0000'),'&#xa;')"/>
        </xsl:for-each> 
        <xsl-text>#;#</xsl-text>  
                       

        <xsl-text>basisPerAtom:</xsl-text>                     
        <xsl:for-each select="$molecule/cml:atomArray/cml:atom">
            <xsl:variable name="outerIndex" select="position()"/>
            <xsl:variable name="elementType" select="@elementType"/>
            <xsl:variable name="id" select="@id"/>              
            <xsl:variable name="basis" select="($fragmentFiles//cml:atom[@elementType=$elementType]/cml:scalar[@dictRef='cc:basis']/text())[1]"/>
            <xsl:variable name="contraction" select="($fragmentFiles//cml:atom[@elementType=$elementType]/cml:scalar[@dictRef='cc:contraction']/text())[1]"/>
            <xsl:value-of select="concat($outerIndex,' ',$elementType,' ',$basis,' ',$contraction,'&#xa;')"/>
        </xsl:for-each>    
                              
    </xsl:template>

    <!-- Charge / Multiplicity section -->
    <xsl:template name="chargemultiplicity">

        <xsl:value-of select="concat('charge:',$charge,'#;#')"/>                            
        <xsl:if test="$multiplicity != ''">
                <xsl:value-of select="concat('multiplicity:',$multiplicity,'#;#')"/>
        </xsl:if>                        
        <xsl:if test="exists($spinpolarization)">
            <xsl:value-of select="concat('spinPolarization:',$spinpolarization,'#;#')"/>                                                           
        </xsl:if>            
    </xsl:template>
    
    <!-- Solvent input section (module) -->
    <xsl:template name="solvatationInfo">
        <xsl:if test="exists($solvation)">
           
            <xsl:variable name="cosmo" select="$solvation/list[@id='cosmo']"/>                                                         
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:ndiv'])">
                <xsl:value-of select="concat('ndiv:',$cosmo/scalar[@dictRef='a:ndiv'],'#;#')"/>
            </xsl:if>
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:nfdiv'])">
                <xsl:value-of select="concat('ndivf:',$cosmo/scalar[@dictRef='a:nfdiv'],'#;#')"/>
            </xsl:if>                
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:rsol'])">    
                <xsl:value-of select="concat('rsolv:',$cosmo/scalar[@dictRef='a:rsol'],'#;#')"/>
                <xsl:value-of select="concat('rsolvUnits:',helper:printUnitSymbol($cosmo/scalar[@dictRef='a:rsol']/@units),'#;#')"/>
            </xsl:if>        
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:rminsolv'])">
                <xsl:value-of select="concat('rminSolv:',$cosmo/scalar[@dictRef='a:rminsolv'],'#;#')"/>
                <xsl:value-of select="concat('rminSolvUnits:',helper:printUnitSymbol($cosmo/scalar[@dictRef='a:rminsolv']/@units),'#;#')"/>
            </xsl:if>     
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:ofac'])">
                <xsl:value-of select="concat('overlapFactorSolv:',$cosmo/scalar[@dictRef='a:ofac'],'#;#')"/>
            </xsl:if>        
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:epsl'])">
                <xsl:value-of select="concat('epsSolv:',$cosmo/scalar[@dictRef='a:epsl'],'#;#')"/>
            </xsl:if>        
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:cosmomethod'])">
                <xsl:value-of select="concat('cosmoMethod:',$cosmo/scalar[@dictRef='a:cosmomethod'],'#;#')"/>
            </xsl:if>      
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:ncix'])">
                <xsl:value-of select="concat('maxIterCharges:',$cosmo/scalar[@dictRef='a:ncix'],'#;#')"/>
            </xsl:if> 
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:ccnv'])">
                <xsl:value-of select="concat('chargeConvergence:',$cosmo/scalar[@dictRef='a:ccnv'],'#;#')"/>
            </xsl:if> 
            <xsl:if test="exists($cosmo/scalar[@dictRef='a:cosmomethod'])">
                <xsl:value-of select="concat('gdsfEmpiricalFactor:',$cosmo/scalar[@dictRef='a:gdsf'],'#;#')"/>
            </xsl:if>     
        </xsl:if>
    </xsl:template>
    

    <!-- Thermochemistry -->
    <xsl:template name="thermochemistry">
        <xsl:variable name="thermochemistry" select="./module[@id='finalization']/propertyList/property[@dictRef='cc:thermochemistry']/module[@cmlx:templateRef='thermochemistry']"/>        
        <xsl:if test="exists($thermochemistry)">
            <xsl:variable name="bonding" select="number((.//cml:module[@id='finalization']/cml:module[@id='otherComponents']/cml:module[@cmlx:templateRef='bonding.energy']/cml:module[@cmlx:templateRef='summary']/cml:scalar[@dictRef='cc:total'])[last()])"/>
                        
            <xsl:variable name="thermochemistryEnergies" select="$thermochemistry/module[@cmlx:templateRef='energies']"/>
            <xsl:variable name="temperature" select="$thermochemistryEnergies/scalar[@dictRef='cc:temp']"/>
            <xsl:value-of select="concat('temperature:',$temperature,'#;#')"/> 
            <xsl:variable name="temperatureUnits" select="$thermochemistryEnergies/scalar[@dictRef='cc:temp']/@units"/>
            <xsl:variable name="entropy" select="$thermochemistryEnergies/list[@cmlx:templateRef='entropy']"/>
            <xsl:variable name="internalEnergy" select="$thermochemistryEnergies/list[@cmlx:templateRef='internalEnergy']"/>
            <xsl:variable name="heat" select="$thermochemistryEnergies/list[@cmlx:templateRef='heat']"/>
                
            <xsl:value-of select="concat('entropy:',$entropy/scalar[@dictRef='cc:total'],'#;#')"/>
            <xsl:value-of select="concat('entropyTrans:',$entropy/scalar[@dictRef='cc:transl'],'#;#')"/>
            <xsl:value-of select="concat('entropyRot:',$entropy/scalar[@dictRef='cc:rotat'],'#;#')"/>
            <xsl:value-of select="concat('entropyVib:',$entropy/scalar[@dictRef='cc:vibrat'],'#;#')"/>
            
            <xsl:value-of select="concat('entropyUnits',$entropy/scalar[@dictRef='cc:total']/@units,'#;#')"/>
            
            <xsl:value-of select="concat('thermalEnergyCorr:',$internalEnergy/scalar[@dictRef='cc:total'],'#;#')"/>
            <xsl:value-of select="concat('thermalEnergyTrans:',$internalEnergy/scalar[@dictRef='cc:transl'],'#;#')"/>
            <xsl:value-of select="concat('thermalEnergyRot:',$internalEnergy/scalar[@dictRef='cc:rotat'],'#;#')"/>
            <xsl:value-of select="concat('thermalEnergyVib:',$internalEnergy/scalar[@dictRef='cc:vibrat'],'#;#')"/>
            
            <xsl:value-of select="concat('thermalEnergyUnits',$internalEnergy/scalar[@dictRef='cc:total']/@units,'#;#')"/>
            <xsl:value-of select="concat('thermalContributionUnits',$internalEnergy/scalar[@dictRef='cc:transl']/@units,'#;#')"/>
            
            <xsl:value-of select="concat('totalHeat:',$heat/scalar[@dictRef='cc:total'],'#;#')"/>
            <xsl:value-of select="concat('translationalHeat:',$heat/scalar[@dictRef='cc:transl'],'#;#')"/>
            <xsl:value-of select="concat('rotationalHeat:',$heat/scalar[@dictRef='cc:rotat'],'#;#')"/>
            <xsl:value-of select="concat('vibrationalHeat:',$heat/scalar[@dictRef='cc:vibrat'],'#;#')"/>
           
            <xsl:variable name="internalFinal" select="number((./cml:list[@cmlx:templateRef='internalEnergy']/cml:scalar[@dictRef='cc:total'])[last()])" />
            <xsl:variable name="entropyFinal" select="number((./cml:list[@cmlx:templateRef='entropy']/cml:scalar[@dictRef='cc:total'])[last()])"/>
            <xsl:if test="exists($bonding) and exists($internalFinal) and exists($entropyFinal)">
                 <xsl:variable name="thermalEnergy" select="(($bonding * 96.48531 * 1000 +  $internalFinal * 4184 ))"/>
                 <xsl:variable name="enthalpy" select="(($bonding * 96.48531 * 1000 +  $internalFinal * 4184 ) + (8.31441 * number($temperature)))"/>
                 <xsl:variable name="gibbsFreeEnergy" select="(($bonding * 96.48531 * 1000 +  $internalFinal * 4184 ) + (8.31441 * number($temperature)) - (number($temperature) * $entropyFinal * 4.184)) div 1000"/>     
                 <xsl:value-of select="concat('thermalEnergy:',round-half-to-even($thermalEnergy div 4.184 ,1),'#;#')"/>
                 <xsl:value-of select="concat('thermalEnergyUnits:','nonsi2.kcal.mol-1','#;#')"/>
                 <xsl:value-of select="concat('enthalpy:',round-half-to-even($enthalpy div 4.184 ,1),'#;#')"/>
                 <xsl:value-of select="concat('enthalpyUnits:','nonsi2.kcal.mol-1','#;#')"/>  
                 <xsl:value-of select="concat('gibbsFreeEnergy:',round-half-to-even($gibbsFreeEnergy div 4.184 ,1),'#;#')"/>
                 <xsl:value-of select="concat('gibbsFreeEnergyUnits:','nonsi2.kcal.mol-1','#;#')"/>
            </xsl:if> 
        <!--            rest of parameters-->
        <xsl:variable name="symmnumb" select="$thermochemistry/cml:scalar[@dictRef='cc:symmnumber']"/>
        <xsl:variable name="mominertia" select="$thermochemistry/cml:array[@dictRef='cc:moi']"/>
        <xsl:variable name="mominertiaUnits" select="'atomic_units'"/>
         
        <xsl:variable name="atomMasses" select="//cml:module[@cmlx:templateRef='adf.frequencyanalysis']/cml:module[@cmlx:templateRef='masses']/cml:array[@dictRef='cc:atomicmass']"/>
        <xsl:variable name="atomMassesProc" select="for $i in tokenize($atomMasses,'\s+') return number($i)"/>
        <xsl:variable name="molmass" select="sum($atomMassesProc)"/>
        <xsl:variable name="molmassUnits" select="'nonsi_amu'"/>
        <xsl:value-of select="concat('molmass:',$molmass,'#;#')"/>
        <xsl:value-of select="concat('molmassUnits:',$molmassUnits,'#;#')"/>    
            
        </xsl:if>
    </xsl:template>
    
    <!-- Intensities -->
    <xsl:template name="frequencyIntensities">
        <xsl:variable name="intensities" select="adf:getIntensities(
                                                    .//module[@id='finalization']//property[@dictRef='cc:intensities']/module[@cmlx:templateRef='intensities'],
                                                    .//module[@id='finalization']//property/module[@cmlx:templateRef='scanfreq']
                                                    )"/>
        <xsl:if test="exists($intensities)">
            <xsl:variable name="frequency" select="tokenize($intensities/array[@dictRef='cc:frequency'],'\s+')"/>
            <xsl:variable name="dipole" select="tokenize($intensities/array[@dictRef='cc:dipole'],'\s+')"/>
            <xsl:variable name="absorption" select="tokenize($intensities/array[@dictRef='cc:absortion'],'\s+')"/>
            
<!--            <xsl-text>frequencies:</xsl-text>
            <xsl:for-each select="1 to count($frequency)">
                <xsl:variable name="innerIndex" select="."/>
                <xsl:value-of select="concat($frequency[$innerIndex],' ')"/>                        
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>-->
            <xsl-text>dipoleStrength:</xsl-text>
            <xsl:for-each select="1 to count($frequency)">
                <xsl:variable name="innerIndex" select="."/>
                <xsl:value-of select="concat($dipole[$innerIndex],' ')"/>                   
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>
            <xsl-text>absorptionIntensity:</xsl-text>
            <xsl:for-each select="1 to count($frequency)">
                <xsl:variable name="innerIndex" select="."/>
                <xsl:value-of select="concat($absorption[$innerIndex],' ')"/>                            
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>
            
            <xsl:value-of select="concat('frequencyUnits:',helper:printUnitSymbol($intensities/array[@dictRef='cc:frequency']/@units),'#;#')"/>
            <xsl:value-of select="concat('dipoleStrengthUnits:',helper:printUnitSymbol($intensities/array[@dictRef='cc:dipole']/@units),'#;#')"/>
            <xsl:value-of select="concat('absorptionIntensityUnits:',helper:printUnitSymbol($intensities/array[@dictRef='cc:absortion']/@units),'#;#')"/>
                      
        </xsl:if>
    </xsl:template>
    
    <xsl:template name="frequencyDisplacements">
        <xsl:variable name="vibrDispl" select="//module[@id='finalization']//property[@dictRef='cc:frequencies']/module[@cmlx:templateRef='vibrations']/cml:array[@dictRef='cc:displacement']"/>       
        <xsl:value-of select="concat('vibrdispl: ',$vibrDispl,'#;#')"/>
        <xsl:text> vibrdisplUnits: nonsi_angstrom #;# </xsl:text>
    </xsl:template>
    
<!--    <!-\- NMR -\->
    <xsl:template name="nmr">
        <xsl:variable name="nmr" select=".//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='nmr']"/>
        <xsl:if test="exists($nmr)">
            <xsl:variable name="elementType" select="tokenize($nmr//array[@dictRef='cc:elementType'], '\s+')"/>
            <xsl:variable name="nucleus" select="tokenize($nmr//array[@dictRef='a:nucleus'], '\s+')"/>
            <xsl:variable name="paramagnetic" select="tokenize($nmr//array[@dictRef='a:paramagneticShielding'], '\s+')"/>
            <xsl:variable name="diamagnetic" select="tokenize($nmr//array[@dictRef='a:diamagneticShielding'], '\s+')"/>
            <xsl:variable name="spinOrbit" select="tokenize($nmr//array[@dictRef='a:spinorbitShielding'], '\s+')"/>
            <xsl:variable name="total" select="tokenize($nmr//array[@dictRef='a:total'], '\s+')"/>
            
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#nmr-{generate-id($nmr)}" style="cursor: pointer; cursor: hand;">
                    <h4 class="panel-title">
                        NMR Shielding Tensors                           
                    </h4>
                </div>
                <div id="nmr-{generate-id($nmr)}" class="panel-collapse collapse">
                    <div class="panel-body">    
                        <div class="row bottom-buffer">
                            <div class="col-md-12 col-sm-12">                           
                                <table id="nmrT-{generate-id($nmr)}">
                                    <thead>
                                        <tr>
                                            <th class="right">Atom</th>   
                                            <th class="right">Paramagnetic <xsl:value-of select="concat('(', helper:printUnitSymbol($nmr//array[@dictRef='a:paramagneticShielding']/@units),')')"/></th>
                                            <th class="right">Diamagnetic <xsl:value-of select="concat('(', helper:printUnitSymbol($nmr//array[@dictRef='a:diamagneticShielding']/@units),')')"/></th>
                                            <th class="right">Spin-orbit <xsl:value-of select="concat('(', helper:printUnitSymbol($nmr//array[@dictRef='a:spinorbitShielding']/@units),')')"/></th>
                                            <th class="right">Total <xsl:value-of select="concat('(', helper:printUnitSymbol($nmr//array[@dictRef='a:total']/@units),')')"/></th>
                                        </tr>                                        
                                    </thead>   
                                    <tbody>
                                        <xsl:for-each select="1 to count($elementType)">
                                            <xsl:variable name="outerIndex" select="."/>
                                            <tr>
                                                <td class="right"><xsl:value-of select="concat($elementType[$outerIndex],'(',$nucleus[$outerIndex],')')"/></td>   
                                                <td class="right"><xsl:value-of select="$paramagnetic[$outerIndex]"/></td>
                                                <td class="right"><xsl:value-of select="$diamagnetic[$outerIndex]"/></td>
                                                <td class="right"><xsl:value-of select="$spinOrbit[$outerIndex]"/></td>
                                                <td class="right"><xsl:value-of select="$total[$outerIndex]"/></td>
                                            </tr>
                                        </xsl:for-each>                                        
                                    </tbody>
                                </table>
                                <script type="text/javascript">
                                    $(document).ready(function() {
                                    $('table#nmrT-<xsl:value-of select="generate-id($nmr)"/>').dataTable({
                                            "bFilter": false,
                                            "bPaginate": false,
                                            "bSort": false,
                                            "bInfo": false
                                        });
                                    } );
                                </script>
                            </div>
                        </div>
                    </div>
                </div>
            </div> 
        </xsl:if>
    </xsl:template>-->
 
    <!-- Zero point vibrational energy -->
    <xsl:template name="zeropointenergy">        
        <xsl:variable name="zeroPoint" select=".//module[@id='finalization']/propertyList/property[@dictRef='cc:zeropoint']/scalar"/>               
        <xsl:if test="exists($zeroPoint)">
            <xsl:value-of select="concat('zeroPointEnergy:',$zeroPoint/text(),'#;#')"/>
            <xsl:value-of select="concat('zeroPointEnergyUnits:',helper:printUnitSymbol($zeroPoint/@units),'#;#')"/>
            
        </xsl:if>
    </xsl:template>

    <!-- Bonding energy -->
    <xsl:template name="bondingEnergy">
        <xsl:variable name="summary" select="(.//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='bonding.energy']/module[@cmlx:templateRef='summary'])[last()]"/>
        <xsl:if test="exists($summary)">
            <xsl-text>energyDecomposition: electrostaticEnergy kineticEnergy coulombEnergy XCEnergy solvationEnergy dispersionEnergy bondingEnergy &#xa;</xsl-text>
            <xsl:value-of select="$summary/scalar[@dictRef='cc:eener']"/><xsl-text> </xsl-text>
            <xsl:value-of select="$summary/scalar[@dictRef='cc:kinener']"/><xsl-text> </xsl-text>
            <xsl:value-of select="$summary/scalar[@dictRef='cc:coulombener']"/><xsl-text> </xsl-text>
            <xsl:value-of select="$summary/scalar[@dictRef='cc:xcener']"/><xsl-text> </xsl-text>
            <xsl:choose>
                <xsl:when test="exists($summary/scalar[@dictRef='cc:solvener'])">
                    <xsl:value-of select="$summary/scalar[@dictRef='cc:solvener']"/><xsl-text> </xsl-text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:text>None </xsl:text>
                </xsl:otherwise>
            </xsl:choose>
            <xsl:choose>
                <xsl:when test="exists($summary/scalar[@dictRef='cc:dispener'])">
                    <xsl:value-of select="$summary/scalar[@dictRef='cc:dispener']"/><xsl-text> </xsl-text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:text>None </xsl:text>
                </xsl:otherwise>
            </xsl:choose>
            <xsl:value-of select="$summary/scalar[@dictRef='cc:total']"/><xsl-text> </xsl-text>
            <xsl-text>#;#</xsl-text>
            <xsl:value-of select="concat('electrostaticEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:eener']/@units),'#;#')"/>
            <xsl:value-of select="concat('kineticEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:kinener']/@units),'#;#')"/>
            <xsl:value-of select="concat('coulombEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:coulombener']/@units),'#;#')"/>
            <xsl:value-of select="concat('XCEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:xcener']/@units),'#;#')"/>
            <xsl:value-of select="concat('electronicEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:total']/@units),'#;#')"/>
            <xsl:value-of select="concat('electronicEnergy:',$summary/scalar[@dictRef='cc:total'],'#;#')"/>
            <xsl:if test="exists($summary/scalar[@dictRef='cc:solvener'])">
                <xsl:value-of select="concat('solvationEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:solvener']/@units),'#;#')"/>
            </xsl:if>
            <xsl:if test="exists($summary/scalar[@dictRef='cc:dispener'])">
                <xsl:value-of select="concat('dispersionEnergyUnits:',helper:printUnitSymbol($summary/scalar[@dictRef='cc:dispener']/@units),'#;#')"/>
            </xsl:if>
                          
        </xsl:if>        
    </xsl:template>

<!--    <!-\- Molecular Orbitals -\->     
    <xsl:template name="molecularOrbitals">
        <xsl:variable name="molecularOrbitals" select="(.//module[@id='finalization']//module[@cmlx:templateRef='sfo.population']/module[@cmlx:templateRef='molecular.orbitals'])[last()]"/>
        <xsl:if test="exists($molecularOrbitals)">
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#molecularOrbitals-{generate-id($molecularOrbitals)}" style="cursor: pointer; cursor: hand;">
                    <h4 class="panel-title">
                        MOs / SFO gross populations                           
                    </h4>
                </div>
                <div id="molecularOrbitals-{generate-id($molecularOrbitals)}" class="panel-collapse collapse">
                    <div class="panel-body">    
                        <div class="row bottom-buffer">
                            <div class="col-lg-8 col-md-8 col-sm-12">
                                <input type="radio" value="all" name="orbitalrangetype" onclick="javascript:filterTable(0)" checked="true">All</input>
                                <input type="radio" value="homolumo" name="orbitalrangetype" onclick="javascript:filterHomoLumo()">Homo/Lumo range:</input>
                                <input type="text"  value="10"       id="orbitalrange" oninput="javascript:filterHomoLumo()" />                                
                                <table id="molecularOrbitalsT-{generate-id($molecularOrbitals)}"></table>
                                <script type="text/javascript">
                                    var molecularOrbitalsData = [  
                                    <xsl:for-each select="$molecularOrbitals/list[@cmlx:templateRef='mo']">
                                        <xsl:variable name="outerIndex" select="position()"/> 
                                        <xsl:variable name="energy" select="./scalar[@dictRef='cc:mo.energy']"/>
                                        <xsl:variable name="occupation" select="./scalar[@dictRef='cc:mo.occupation']"/>
                                        <xsl:variable name="number" select="./scalar[@dictRef='cc:mo.number']"/>
                                        <xsl:variable name="label" select="./scalar[@dictRef='cc:mo.label']"/>
                                        <xsl:variable name="percent" select="tokenize(./array[@dictRef='cc:percent'],'\s+')"/>
                                        <xsl:variable name="sfo1" select="tokenize(./array[@dictRef='cc:sfo1'],'\s+')"/>
                                        <xsl:variable name="sfo2" select="tokenize(./array[@dictRef='cc:sfo2'],'\s+')"/>
                                        <xsl:variable name="partialEnergy" select="tokenize(./array[@dictRef='cc:energy'],'\s+')"/>
                                        <xsl:variable name="partialOccupation" select="tokenize(./array[@dictRef='cc:occupation'],'\s+')"/>
                                        <xsl:variable name="fragment1" select="tokenize(./array[@dictRef='cc:fragment1'],'\s+')"/>
                                        <xsl:variable name="fragment2" select="tokenize(./array[@dictRef='cc:fragment2'],'\s+')"/>                                        
                                        <xsl:for-each select="1 to count($partialEnergy)">
                                            <xsl:variable name="innerIndex" select="."/>
                                            <xsl:choose>
                                                <xsl:when test="$innerIndex = 1">[ <xsl:value-of select="$energy"/>, <xsl:value-of select="$occupation"/>,<xsl:value-of select="$number"/>, "<xsl:value-of select="$label"/>", <xsl:value-of select="$percent[$innerIndex]"/>, <xsl:value-of select="$sfo1[$innerIndex]"/>, "<xsl:value-of select="$sfo2[$innerIndex]"/>", <xsl:value-of select="$partialEnergy[$innerIndex]"/>,<xsl:value-of select="$partialOccupation[$innerIndex]"/>,<xsl:value-of select="$fragment1[$innerIndex]"/>, "<xsl:value-of select="$fragment2[$innerIndex]"/>"]<xsl:if test="($innerIndex &lt; count($partialEnergy)) or ($outerIndex &lt; count($molecularOrbitals/list[@cmlx:templateRef='mo']))"><xsl:text>,</xsl:text></xsl:if></xsl:when>
                                                <xsl:otherwise>["", "", "", "", <xsl:value-of select="$percent[$innerIndex]"/>, <xsl:value-of select="$sfo1[$innerIndex]"/>, "<xsl:value-of select="$sfo2[$innerIndex]"/>", <xsl:value-of select="$partialEnergy[$innerIndex]"/>,<xsl:value-of select="$partialOccupation[$innerIndex]"/>,<xsl:value-of select="$fragment1[$innerIndex]"/>, "<xsl:value-of select="$fragment2[$innerIndex]"/>"]<xsl:if test="($innerIndex &lt; count($partialEnergy)) or ($outerIndex &lt; count($molecularOrbitals/list[@cmlx:templateRef='mo']))"><xsl:text>,</xsl:text></xsl:if></xsl:otherwise>
                                            </xsl:choose>                                                                    
                                        </xsl:for-each>                                                   
                                    </xsl:for-each>
                                    ]; 
                                    
                                    $(document).ready(function() {                        
                                        $('table#molecularOrbitalsT-<xsl:value-of select="generate-id($molecularOrbitals)"/>').dataTable( {
                                            "aaData": molecularOrbitalsData,
                                            "aoColumns": [
                                                { "sTitle": "E(eV)", "sClass" : "right" },
                                                { "sTitle": "Occ"},
                                                { "sTitle": "MO"},
                                                { "sTitle": ""  },
                                                { "sTitle": "%"},
                                                { "sTitle": "SFO" },
                                                { "sTitle": "(first member)" },
                                                { "sTitle": "E(eV)", "sClass" : "right"  },
                                                { "sTitle": "Occ", "sClass" : "right"  },
                                                { "sTitle": "Fragment" },
                                                { "sTitle": "" }
                                            ],
                                            "bFilter": false,
                                            "bPaginate": true,
                                            "bSort": false,
                                            "bInfo": true
                                        } );   
                                    } );                                       
                                </script>                
                                <script type="text/javascript">                                                               
                                    var homoLumoFrontier = null;
                                    var molecularOrbitalsDataShort = null;
                                    
                                    function filterHomoLumo(){
                                        $("input[name='orbitalrangetype'][value='homolumo']").prop('checked',true);
                                        var range = $('input#orbitalrange').val();
                                        if(!$.isNumeric(range)){               
                                            range = 10;
                                        }
                                        filterTable(range);               
                                    }
                                    
                                    function findHomoRange(startIndex, range){
                                        energyChangesCounter = 0;
                                        for(inx = startIndex - 1; inx > 0 ; inx-\- ){
                                            if(molecularOrbitalsData[inx][1].length != 0)                                          
                                                energyChangesCounter++;                                                                                                                    
                                            if(energyChangesCounter == range)
                                                return homoLumoFrontier - inx ;                                            
                                        }
                                        return 0;                                                                        
                                    }
                                    
                                    function findLumoRange(startIndex, range){
                                        energyChangesCounter = 0;    
                                        for(inx = startIndex + 1; inx &#60; molecularOrbitalsData.length ; inx++ ){
                                            if(molecularOrbitalsData[inx][1].length != 0)             
                                                energyChangesCounter++;                                                                                            
                                            if(energyChangesCounter == range)
                                                return (inx -1) - homoLumoFrontier;                                            
                                        }
                                        return molecularOrbitalsData.length - 1;
                                    }                                    
                                    
                                    function filterTable(range){                                                                                                            
                                        findHomoLumoFrontier();
                                        if(homoLumoFrontier == -1 || range == 0){
                                            molecularOrbitalsDataShort = molecularOrbitalsData;
                                        }else{                                        
                                            var lumoRange = findLumoRange(homoLumoFrontier, range);
                                            var homoRange = findHomoRange(homoLumoFrontier, range);
                                        
                                            var homoOrbitals = molecularOrbitalsData.slice(homoLumoFrontier - Number(homoRange) , homoLumoFrontier );
                                            var lumoOrbitals = molecularOrbitalsData.slice(homoLumoFrontier , homoLumoFrontier + Number(lumoRange) + 1);
                                            molecularOrbitalsDataShort = homoOrbitals.concat(lumoOrbitals);
                                        }                                    
                                        var table = $('#molecularOrbitalsT-<xsl:value-of select="generate-id($molecularOrbitals)"/>').DataTable();
                                        table.clear();
                                        table.rows.add(molecularOrbitalsDataShort);
                                        table.draw();                        
                                    }
                                    
                                    function findHomoLumoFrontier(){
                                        if(homoLumoFrontier == null){
                                            homoLumoFrontier = -1;
                                            for(inx = 0; inx <xsl:text disable-output-escaping="yes">&lt;</xsl:text> molecularOrbitalsData.length; inx++){
                                                if( (molecularOrbitalsData[inx][1] == 0) <xsl:text disable-output-escaping="yes">&amp;&amp;</xsl:text> (molecularOrbitalsData[inx][1].length != 0 )){
                                                    homoLumoFrontier = inx;
                                                    break;
                                                }
                                            }                              
                                        }
                                        return homoLumoFrontier;                          
                                    }
                                </script>                                                       
                            </div>                
                        </div>
                    </div>
                </div>
            </div>
        </xsl:if>
    </xsl:template>-->

    <!-- Squared S -->
    <xsl:template name="s2">
        <xsl:variable name="s2" select="(.//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='s2'])"/>
        <xsl:variable name="s2exact" select="$s2/scalar[@dictRef='cc:s2']"/>
        <xsl:variable name="s2expected" select="$s2/scalar[@dictRef='cc:s2expected']"/>
        <xsl:if test="exists($s2)">
            <xsl:value-of select="concat('s2Exact:',$s2exact/text(),'#;#')"></xsl:value-of>
            <xsl:value-of select="concat('s2Expected:',$s2expected/text(),'#;#')"></xsl:value-of>     
        </xsl:if> 
    </xsl:template>
    
    <!-- Fit Test -->
    <xsl:template name="fitTest">
        <xsl:variable name="fitTest" select="(.//module[@id='finalization']/module[@id='otherComponents']/module[@cmlx:templateRef='fit.test'])[last()]"/>        
        <xsl:if test="exists($fitTest)">
            <xsl:value-of select="concat('sumFragments:',$fitTest/scalar[@dictRef='cc:sumfragments']/text(),'#;#')"/>
            <xsl:value-of select="concat('orthogonalizedFragments:',$fitTest/scalar[@dictRef='cc:ortho']/text(),'#;#')"/>
            <xsl:value-of select="concat('fitSCF:',$fitTest/scalar[@dictRef='cc:fitscf']/text(),'#;#')"/>
        </xsl:if> 
    </xsl:template>
    
    <!-- Mulliken charges and orbital occupation -->
    <xsl:template name="mullikenCharges">
        <xsl:variable name="charges" select="(.//module[@id='finalization']//module[@cmlx:templateRef='mulliken']/module[@cmlx:templateRef='charges'])[last()]"/>
        <xsl:if test="exists($charges)">
            <xsl:variable name="serial" select="tokenize($charges/list/*[@dictRef='cc:serial'],'\s+')"/>
            <xsl:variable name="elementType" select="tokenize($charges/list/*[@dictRef='cc:elementType'],'\s+')"/>
            <xsl:variable name="charge" select="tokenize($charges/list/*[@dictRef='x:charge'],'\s+')"/>
            <xsl:variable name="spinDensity" select="tokenize($charges/list/*[@dictRef='a:spinDensity'],'\s+')"/>
            <xsl:variable name="spin" select="tokenize($charges/list/*[@dictRef='a:spin'],'\s+')"/>
            <xsl:variable name="orbitalS" select="tokenize($charges/list/*[@dictRef='a:orbitalS'],'\s+')"/>
            <xsl:variable name="orbitalP" select="tokenize($charges/list/*[@dictRef='a:orbitalP'],'\s+')"/>
            <xsl:variable name="orbitalD" select="tokenize($charges/list/*[@dictRef='a:orbitalD'],'\s+')"/>
            <xsl:variable name="orbitalF" select="tokenize($charges/list/*[@dictRef='a:orbitalF'],'\s+')"/>            
            
            <xsl:choose> 
                <xsl:when test="not(exists($spin))">
                    <xsl-text>mullikenCharges: index atom charge orbS orbP orbD orbF &#xa; </xsl-text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl-text>mullikenCharges: index atom charge spinDensity spin_A orbS_A orbP_A orbD_A orbF_A spin_B orbS_B orbP_B orbD_B orbF_B &#xa; </xsl-text>               
                </xsl:otherwise>
            </xsl:choose>
            
            <xsl:for-each select="1 to count($serial)">
                <xsl:variable name="innerIndex" select="."/>                                        
                <xsl:choose>
                    <xsl:when test="not(exists($spin))">
                        <xsl:value-of select="$serial[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$elementType[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$charge[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$orbitalS[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$orbitalP[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$orbitalD[$innerIndex]"/> <xsl:text> </xsl:text>
                        <xsl:value-of select="$orbitalF[$innerIndex]"/> <xsl:text> &#xa; </xsl:text>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:value-of select="$serial[$innerIndex]"/>
                        <xsl:value-of select="$elementType[$innerIndex]"/>
                        <xsl:value-of select="$charge[$innerIndex]"/>
                        <xsl:value-of select="$spinDensity[$innerIndex]"/>
                        <xsl:value-of select="$spin[$innerIndex * 2 - 1]"/>
                        <xsl:value-of select="$orbitalS[$innerIndex * 2 - 1]"/>
                        <xsl:value-of select="$orbitalP[$innerIndex * 2 - 1]"/>
                        <xsl:value-of select="$orbitalD[$innerIndex * 2 - 1]"/>
                        <xsl:value-of select="$orbitalF[$innerIndex * 2 - 1]"/>"
                        <xsl:value-of select="$spin[$innerIndex * 2]"/>
                        <xsl:value-of select="$orbitalS[$innerIndex * 2]"/>
                        <xsl:value-of select="$orbitalP[$innerIndex * 2]"/>
                        <xsl:value-of select="$orbitalD[$innerIndex * 2]"/>
                        <xsl:value-of select="$orbitalF[$innerIndex * 2]"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:for-each>      
            <xsl-text>#;#</xsl-text>
                                  
        </xsl:if>
    </xsl:template>
    
    <!-- Frequencies section -->    
    <xsl:template name="frequencies">
        <xsl:variable name="frequencies" select="adf:getFrequencies(
                                                        .//module[@id='finalization']/propertyList/property[@dictRef='cc:frequencies']/module[@dictRef='cc:vibrations'],
                                                        .//module[@id='finalization']//property/module[@cmlx:templateRef='scanfreq']
                                                        )"/>       
        <xsl:if test="exists($frequencies)">
            <xsl:variable name="frequency" select="tokenize($frequencies/array[@dictRef='cc:frequency'],'\s+')"/>
            <xsl-text>frequencies:</xsl-text>
            <xsl:for-each select="1 to count($frequency)">
                <xsl:variable name="innerIndex" select="."/>
                <xsl:value-of select="concat($frequency[$innerIndex],' ')"/>                        
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>
            <xsl:text>freqUnits:nonsi_cm-1#;#</xsl:text>
        </xsl:if>
    </xsl:template>
       
    <!-- Multipole derived atomic charges -->
    <xsl:template name="mdc">
        <xsl:variable name="atomicCharges" select="(.//module[@id='finalization']//module[@cmlx:templateRef='atomic.charges']/list[@cmlx:templateRef='multipole'])[last()]"/>
        <xsl:variable name="atomicChargesSpin" select="(.//module[@id='finalization']//module[@cmlx:templateRef='atomic.charges.spin']/list[@cmlx:templateRef='spin'])[last()]"/>
        <xsl:variable name="spinDensity" select="(.//module[@id='finalization']//module[@cmlx:templateRef='spin.density']/list[@cmlx:templateRef='spinDensity'])[last()]"/>
        <xsl:variable name="uid" select="concat(generate-id($atomicCharges),generate-id($atomicChargesSpin),generate-id($spinDensity))"/>
        <xsl:if test="exists($atomicCharges) or exists($atomicChargesSpin) or exists($spinDensity)">
                            
             <xsl:call-template name="atomicCharges">
                 <xsl:with-param name="atomicCharges" select="$atomicCharges"/>
             </xsl:call-template>
             <xsl:call-template name="atomicChargesSpin">
                <xsl:with-param name="atomicChargesSpin" select="$atomicChargesSpin"/>
             </xsl:call-template>                              
             <xsl:call-template name="spinDensity">
                <xsl:with-param name="spinDensity" select="$spinDensity"/>
             </xsl:call-template>
                      
        </xsl:if>
    </xsl:template>
    <!-- Multipole derived atomic charges : atomic charges -->
    <xsl:template name="atomicCharges">
        <xsl:param name="atomicCharges"/>
        <xsl:if test="exists($atomicCharges)">
            <xsl:variable name="serial" select="tokenize($atomicCharges/array[@dictRef='cc:serial'], '\s+')"/>
            <xsl:variable name="elementType" select="tokenize($atomicCharges/array[@dictRef='cc:elementType'], '\s+')"/>
            <xsl:variable name="mdcm" select="tokenize($atomicCharges/array[@dictRef='a:mdcm'], '\s+')"/>
            <xsl:variable name="mdcd" select="tokenize($atomicCharges/array[@dictRef='a:mdcd'], '\s+')"/>
            <xsl:variable name="mdcq" select="tokenize($atomicCharges/array[@dictRef='a:mdcq'], '\s+')"/>
            <xsl-text>multipoleDerivedCharges: index atom mdcM mdcD mdcQ &#xa;</xsl-text>
            <xsl:for-each select="1 to count($serial)">
                <xsl:variable name="innerIndex" select="."/>                                                          
                <xsl:value-of select="concat($serial[$innerIndex],' ')"/>
                <xsl:value-of select="concat($elementType[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcm[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcd[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcq[$innerIndex],'&#xa;')"/>   
            </xsl:for-each> 
            <xsl-text>#;#</xsl-text>       
        </xsl:if>
    </xsl:template>
    <!-- Multipole derived atomic charges : atomic charges spin -->
    <xsl:template name="atomicChargesSpin">
        <xsl:param name="atomicChargesSpin"/>            
        <xsl:if test="exists($atomicChargesSpin)">
            <xsl:variable name="serial" select="tokenize($atomicChargesSpin/array[@dictRef='cc:serial'], '\s+')"/>
            <xsl:variable name="elementType" select="tokenize($atomicChargesSpin/array[@dictRef='cc:elementType'], '\s+')"/>
            <xsl:variable name="mdcm" select="tokenize($atomicChargesSpin/array[@dictRef='a:mdcm'], '\s+')"/>
            <xsl:variable name="mdcd" select="tokenize($atomicChargesSpin/array[@dictRef='a:mdcd'], '\s+')"/>
            <xsl:variable name="mdcq" select="tokenize($atomicChargesSpin/array[@dictRef='a:mdcq'], '\s+')"/>
            <xsl-text>multipoleDerivedChargesSpinAB: index atom mdcM mdcD mdcQ &#xa;</xsl-text>
            <xsl:for-each select="1 to count($serial)">
                <xsl:variable name="innerIndex" select="."/>                                                          
                <xsl:value-of select="concat($serial[$innerIndex],' ')"/>
                <xsl:value-of select="concat($elementType[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcm[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcd[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcq[$innerIndex],'&#xa;')"/>              
            </xsl:for-each>      
            <xsl-text>#;#</xsl-text> 
        </xsl:if>
    </xsl:template>
    <!-- Multipole derived atomic charges : spin density -->
    <xsl:template name="spinDensity">
        <xsl:param name="spinDensity"/>        
        <xsl:if test="exists($spinDensity)">
            <xsl:variable name="serial" select="tokenize($spinDensity/array[@dictRef='cc:serial'], '\s+')"/>
            <xsl:variable name="elementType" select="tokenize($spinDensity/array[@dictRef='cc:elementType'], '\s+')"/>
            <xsl:variable name="mdcm" select="tokenize($spinDensity/array[@dictRef='a:mdcm'], '\s+')"/>
            <xsl:variable name="mdcd" select="tokenize($spinDensity/array[@dictRef='a:mdcd'], '\s+')"/>
            <xsl:variable name="mdcq" select="tokenize($spinDensity/array[@dictRef='a:mdcq'], '\s+')"/>
            <xsl-text>spinDensity: index atom mdcM mdcD mdcQ &#xa;</xsl-text>
        
            <xsl:for-each select="1 to count($serial)">
                <xsl:variable name="innerIndex" select="."/>                                                          
                <xsl:value-of select="concat($serial[$innerIndex],' ')"/>
                <xsl:value-of select="concat($elementType[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcm[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcd[$innerIndex],' ')"/>
                <xsl:value-of select="concat($mdcq[$innerIndex],'&#xa;')"/>  
            </xsl:for-each>    
            <xsl-text>#;#</xsl-text> 
        </xsl:if>
    </xsl:template>
    
    <!-- Orbital energies -->
    <xsl:template name="orbitalEnergy">
        <xsl:variable name="orbitalEnergies" select="
            if(exists(.//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies']/list[@cmlx:templateRef='energies'])) then
                .//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies']/list[@cmlx:templateRef='energies']
            else
                .//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies.zora']/list[@cmlx:templateRef='energies']
            "/>
        
        <xsl:variable name="orbitalEnergiesSpin" select="
            if(exists(.//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies.spin']/list[@cmlx:templateRef='energies'])) then 
                .//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies.spin']/list[@cmlx:templateRef='energies']
            else
                .//module[@id='finalization']//module[@cmlx:templateRef='orbital.energies.spin.zora']/list[@cmlx:templateRef='energies']
            "/>
        <xsl:variable name="uid" select="concat(generate-id($orbitalEnergies),generate-id($orbitalEnergiesSpin))"/>
        <xsl:if test="exists($orbitalEnergies) or exists($orbitalEnergiesSpin)">

          <xsl:call-template name="orbitalEnergies">
              <xsl:with-param name="orbitalEnergies" select="$orbitalEnergies"/>
          </xsl:call-template>
          <xsl:call-template name="orbitalEnergiesSpin">
              <xsl:with-param name="orbitalEnergiesSpin" select="$orbitalEnergiesSpin"/>
          </xsl:call-template>                              
                           
        </xsl:if>     
    </xsl:template>
    <xsl:template name="orbitalEnergies">
        <xsl:param name="orbitalEnergies"/>
        <xsl:if test="exists($orbitalEnergies)">
            <xsl:variable name="irrep" select="tokenize($orbitalEnergies/array[@dictRef='cc:irrep'], '\s+')"/>
            <xsl:variable name="serial" select="tokenize($orbitalEnergies/array[@dictRef='cc:serial'], '\s+')"/>
            <xsl:variable name="occup" select="tokenize($orbitalEnergies/array[@dictRef='cc:occup'], '\s+')"/>
            <xsl:variable name="energy" select="tokenize($orbitalEnergies/array[@dictRef='cc:energy'], '\s+')"/> 
            <xsl-text>orbitalEnergies: irrep index occupation energy &#xa;</xsl-text>                                                               
            <xsl:for-each select="1 to count($irrep)">
                    <xsl:variable name="innerIndex" select="."/>                                                          
                     <xsl:value-of select="concat($irrep[$innerIndex],' ')"/>
                     <xsl:value-of select="concat($serial[$innerIndex],' ')"/>
                     <xsl:value-of select="concat($occup[$innerIndex],' ')"/>
                     <xsl:value-of select="concat($energy[$innerIndex],'&#xa;')"/>
           
            </xsl:for-each>                                     
            <xsl-text>#;#</xsl-text>                       
        </xsl:if>                
    </xsl:template>
    <xsl:template name="orbitalEnergiesSpin">
        <xsl:param name="orbitalEnergiesSpin"/>     
        <xsl:if test="exists($orbitalEnergiesSpin)">
              
            <xsl:variable name="irrep" select="tokenize($orbitalEnergiesSpin/array[@dictRef='cc:irrep'], '\s+')"/>
            <xsl:variable name="serial" select="tokenize($orbitalEnergiesSpin/array[@dictRef='cc:serial'], '\s+')"/>
            <xsl:variable name="spin" select="tokenize($orbitalEnergiesSpin/array[@dictRef='cc:spin'], '\s+')"/>
            <xsl:variable name="occup" select="tokenize($orbitalEnergiesSpin/array[@dictRef='cc:occup'], '\s+')"/>
            <xsl:variable name="energy" select="tokenize($orbitalEnergiesSpin/array[@dictRef='cc:energy'], '\s+')"/>  
            <xsl-text>orbitalEnergiesSpinAB: irrep index occupation spin energy &#xa;</xsl-text>   
            <xsl:for-each select="1 to count($irrep)">
                <xsl:variable name="innerIndex" select="."/>                                                          
                <xsl:value-of select="$irrep[$innerIndex]"/><xsl-text> </xsl-text>
                <xsl:value-of select="$serial[$innerIndex]"/><xsl-text> </xsl-text>
                <xsl:value-of select="$spin[$innerIndex]"/><xsl-text> </xsl-text>
                <xsl:value-of select="$occup[$innerIndex]"/><xsl-text> </xsl-text>
                <xsl:value-of select="$energy[$innerIndex]"/><xsl-text>&#xa;</xsl-text>                 
            </xsl:for-each>                                     
            <xsl-text>#;#</xsl-text>  
          
        </xsl:if> 
    </xsl:template>
    
    <!-- Multipole molecular moments -->
    <xsl:template name="multipoleMoments">
        <xsl:variable name="dipole" select="(.//module[@id='finalization']//module[@cmlx:templateRef='dipole.moment'])[last()]"/>
        <xsl:variable name="quadrupole" select="(.//module[@id='finalization']//module[@cmlx:templateRef='quadrupole.moment']/array[@dictRef='cc:quadrupole'])[last()]"/>        
        <xsl:if test="exists($dipole) or exists($quadrupole)">
            <xsl:if test="exists($dipole)">
                <xsl:value-of select="concat('dipoleMomentVector:',$dipole[1],' ',$dipole[2],' ',$dipole[3],'#;#')"/>
                <xsl:value-of select="concat('dipoleMomentTotal:',$dipole/cml:scalar[@dictRef='x:dipole'],'#;#')"/>
            </xsl:if>
            
            <xsl:if test="exists($quadrupole)">
                <xsl-text>quadrupoleMoment:XX XY ZZ XY XZ YZ &#xa;</xsl-text>
                <xsl:value-of select="$quadrupole"/> 
                <xsl-text>#;#</xsl-text>
            </xsl:if>
        </xsl:if>       
    </xsl:template>
<!--    
    <!-\- Timing -\->
    <xsl:template name="timing">
        <xsl:variable name="cputime" select=".//module[@id='finalization']/propertyList/property[@dictRef='cc:cputime']/scalar"/>
        <xsl:variable name="systemtime" select=".//module[@id='finalization']/propertyList/property[@dictRef='cc:systemtime']/scalar"/>
        <xsl:variable name="elapsedtime" select=".//module[@id='finalization']/propertyList/property[@dictRef='cc:elapsedtime']/scalar"/>                
        <xsl:if test="exists($cputime)">
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#timing-{generate-id($cputime)}" style="cursor: pointer; cursor: hand;">
                    <h4 class="panel-title">
                        Timing                        
                    </h4>
                </div>
                <div id="timing-{generate-id($cputime)}" class="panel-collapse collapse">
                    <div class="panel-body">    
                        <div class="row bottom-buffer">
                            <div class="col-lg-4 col-md-4 col-sm-6"> 
                                <table id="timingT-{generate-id($cputime)}">
                                    <thead>
                                        <tr>
                                            <th>Factor</th>
                                            <th> </th>
                                        </tr> 
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td>Cpu</td>
                                            <td><xsl:value-of select="$cputime"/></td>                        
                                        </tr>                         
                                        <tr>
                                            <td>System</td>
                                            <td><xsl:value-of select="$systemtime"/></td>                        
                                        </tr>                                             
                                        <tr>
                                            <td>Elapsed</td>
                                            <td><xsl:value-of select="$elapsedtime"/></td>                        
                                        </tr>                                             
                                    </tbody>
                                </table>
                                <script type="text/javascript">
                                    $(document).ready(function() {                        
                                    $('table#timingT-<xsl:value-of select="generate-id($cputime)"/>').dataTable( {
                                    "bFilter": false,                                 
                                    "bPaginate": false,                                    
                                    "bSort": false,
                                    "bInfo": false,
                                    "aoColumnDefs" : [                                    
                                        { "sClass": "text-right", "aTargets": [ 1 ] }
                                    ]
                                    } );   
                                    } );                                       
                                </script>  
                            </div>
                        </div>
                    </div>
                </div>
            </div>                                  
        </xsl:if>                
    </xsl:template>-->
    
 <!--   <!-\- Excitations -\->
    <xsl:template name="excitations">
        <xsl:variable name="excitations" select=".//module[@id='finalization']//module[@cmlx:templateRef='excitation.energy']"/>
        <xsl:if test="exists($excitations)">
            <xsl:variable name="excitationEnergies" select="$excitations/module[@cmlx:templateRef='excitationEnergies']"/>             
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#excitations-{generate-id($excitations)}" style="cursor: pointer; cursor: hand;">                    
                    <h4 class="panel-title">
                        Final Excitation Energies                        
                    </h4>
                </div>
                <div id="excitations-{generate-id($excitations)}" class="panel-collapse collapse">
                    <div class="panel-body">
                        
                        <xsl:variable name="bandwidthev" select="0.15"/>
                        <xsl:variable name="bandwidthcm1" select="1200"/>
                        <xsl:variable name="bandwidthnm" select="20"/>                        
                        <script type="text/javascript">                                    
    
                            var evBandwidth = <xsl:value-of select="$bandwidthev"/>;
                            var cm1Bandwidth   = <xsl:value-of select="$bandwidthcm1"/>;
                            var nmBandwidth = <xsl:value-of select="$bandwidthnm"/>;
                            
                            var evStep = evBandwidth / 30;
                            var cm1Step = cm1Bandwidth / 30;
                            var nmStep = nmBandwidth / 30;
                            
                            var oscStrengths = new Array();                          
                            var energyEv     = new Array();                     
                            var energyCm1    = new Array();
                            var energyNm     = new Array();
                            
                            var energyMaxEv = new Array();                 
                            var energyMinEv = new Array();                
                            var energyMaxCm1 = new Array();          
                            var energyMinCm1 = new Array();
                            var energyMaxNm = new Array();   
                            var energyMinNm = new Array();   
                            
                            function round(num){
                                return Math.round(num * 100) / 100                                    
                            }
                            
                            function excitationParametersChanged(plotId){
                                $("excitation-generateBtn-" + "excitationPlotId").prop('disabled', false);
                                setExcitationGraphLimitValues(plotId);                                                                   
                                generateExcitationChart(plotId);                       
                            }
                            
                            function setExcitationGraphLimitValues(plotId){
                                value = $("input[name='excitationGraphUnits-" + plotId + "']:checked").val();
                                if(value == "ev"){                        		     
                                    $("#maxXEnergy-" + plotId).val(energyMaxEv[plotId]);
                                    $("#minXEnergy-" + plotId).val(energyMinEv[plotId]);      
                                    $("#bandwidth-" + plotId).val(evBandwidth);            
                                }
                                if(value == "cm1"){     
                                    $("#maxXEnergy-" + plotId).val(energyMaxCm1[plotId]);
                                    $("#minXEnergy-" + plotId).val(energyMinCm1[plotId]);
                                    $("#bandwidth-" + plotId).val(cm1Bandwidth);
                                }
                                if(value == "nm"){
                                   $("#maxXEnergy-" + plotId).val(energyMaxNm[plotId]);
                                   $("#minXEnergy-" + plotId).val(energyMinNm[plotId]);
                                   $("#bandwidth-" + plotId).val(nmBandwidth);
                                }
                            }
                            
                            function generateExcitationChart(plotId){
                                $('#excitationEnergyPlotlyContainer-' + plotId).show();
                                var data = fillExcitationChart(plotId);                                    
                                var xLabel = "";                                        
                                unitType = $("input[name='excitationGraphUnits-" + plotId + "']:checked").val();
                                if(unitType == "ev")                        		              
                                    xLabel = 'E (eV)';
                                else if(unitType == "cm1")
                                    xLabel = 'E (cm-1)';
                                else if(unitType == "nm")
                                    xLabel = "nm";
                            
                                var layout = {
                                    height: 550,
                                    xaxis: {
                                        title: xLabel,
                                        showgrid: true,
                                        zeroline: true,
                                        tickmode: 'auto'
                                    },
                                    hovermode: 'closest',
                                    yaxis: {
                                        title: 'Oscillator strengths (a.u.)'
                                    }
                                };                                                
                                Plotly.react('excitationEnergyPlotlyContainer-' + plotId, data, layout,  {responsive: true, showSendToCloud: true});
                            }
                            
                            function fillExcitationChart(plotId){
                                //Clear possible previous values                                        
                                var data = new Array();
                                var dataSum = new Array();
                                
                                //Now we'll generate individual energy and all-energy-sum series and append them in data and dataSum arrays
                                var bandWidth = $("#bandwidth-" + plotId).val();
                                var minX = parseFloat($("#minXEnergy-" + plotId).val());
                                var maxX = parseFloat($("#maxXEnergy-" + plotId).val());
                                
                                unitType = $("input[name='excitationGraphUnits-" + plotId + "']:checked").val();                                       
                                for(var inx = 0; inx <xsl:text disable-output-escaping="yes">&lt;</xsl:text> oscStrengths[plotId].length; inx++){
                                    generateGaussianValuesForPeak(plotId, data,dataSum,  inx+1, unitType, bandWidth, minX, maxX);
                                }                                        
                                
                                var xVal = new Array();
                                var yVal = new Array();
                                var sum = 0;
                                for(var inx = 0; inx <xsl:text disable-output-escaping="yes">&lt;</xsl:text> dataSum.length; inx++){
                                    xVal.push(dataSum[inx][0]);
                                    yVal.push(dataSum[inx][1]);
                                    sum += dataSum[inx][1];
                                }
                                
                                //Check if there is representable data, if not (no osc.strength values) we won't show the graph                                        
                                if(sum == 0){
                                    $('#excitationEnergyPlotlyContainer-' + plotId).hide();
                                    $('#excitationEnergyContainerControls-' + plotId).hide();
                                    return;
                                }
                                
                                //Add sum serie and plot graph
                                var serie = {
                                    x: xVal,
                                    y: yVal,
                                    mode: 'lines',
                                    name: 'Total sum',
                                    marker: {                                    
                                            color: 'rgb(57, 57, 69)'
                                    },
                                    line: {
                                            color: 'rgb(57, 57, 69)'
                                    }                                    
                                };                                                           
                                data.push(serie);
                                
                                // Add bars
                                var energyArray;
                                switch(unitType){
                                    case 'ev':  energyArray = energyEv[plotId];
                                                break;
                                    case 'cm1': energyArray = energyCm1[plotId];
                                                break;
                                    case 'nm':  energyArray = energyNm[plotId];
                                                break;         
                                }      
                                
                                var serieBars =  {
                                    x : energyArray,
                                    y : oscStrengths[plotId],
                                    type: 'bar',
                                    mode: 'markers',                                    
                                    name: 'Intensities',
                                    marker: {                                    
                                             color: 'rgb(224, 97, 97)'
                                            },
                                    line: {
                                             color: 'rgb(224, 97, 97)'
                                          }
                                }
                                data.push(serieBars);
                                return data;
                            }
                            
                            function generateGaussianValuesForPeak(plotId, data, dataSum, energyInx, units, bandWidth, minX, maxX) {
                                var step;
                                var energyValue;                                
                                switch(units){
                                    case 'ev':  step = evStep;
                                                energyValue = energyEv[plotId][energyInx-1];
                                                break;
                                    case 'cm1': step = cm1Step;
                                                energyValue = energyCm1[plotId][energyInx-1];
                                                break;
                                    case 'nm':  step = nmStep;
                                                energyValue = energyNm[plotId][energyInx-1];
                                                break;         
                                }      
                                
                                var inx = 0;
                                var xVal = new Array();
                                var yVal = new Array();
                                for(var x = minX; x <xsl:text disable-output-escaping="yes">&lt;</xsl:text> maxX ; x = x + step ) {  
                                    if(dataSum.length <xsl:text disable-output-escaping="yes">&lt;</xsl:text>= inx)     
                                        dataSum.push([x,oscStrengths[plotId][energyInx-1] * Math.exp( -Math.pow(((2*Math.sqrt(2*Math.log(2)))/bandWidth),2) * Math.pow((x - energyValue),2))]);
                                    else
                                        dataSum[inx][1] += oscStrengths[plotId][energyInx-1] * Math.exp( -Math.pow(((2*Math.sqrt(2*Math.log(2)))/bandWidth),2) * Math.pow((x - energyValue),2))
                                    
                                    xVal.push(x);
                                    yVal.push(oscStrengths[plotId][energyInx-1] * Math.exp( -Math.pow(((2*Math.sqrt(2*Math.log(2)))/bandWidth),2) * Math.pow((x -energyValue),2)));
                                    inx++;
                                }
                                
                                // Individual series are disabled right now
                                /*
                                var serie = {
                                    x: xVal,
                                    y: yVal,
                                    mode: 'lines',
                                    name: 'Energy ' + energyInx,
                                };                                                           
                                data.push(serie); 
                                */
                                
                            }
                        </script>                                                
                        
                        <xsl:for-each select="1 to count($excitationEnergies)">
                            <xsl:variable name="outerIndex" select="."/>
                            <xsl:variable name="energies" select="$excitationEnergies[$outerIndex]/module[@cmlx:templateRef='energies']"/>                            
                            <xsl:variable name="dipole"   select="$excitationEnergies[$outerIndex]/module[@cmlx:templateRef='dipole']"/>
                            <xsl:variable name="sym" select="$excitationEnergies[$outerIndex]/scalar[@dictRef='cc:symm']"></xsl:variable>
                            <!-\- Print header -\->
                           <div class="row">
                               <div class="col-md-12">
                                   <strong>Symmetry <xsl:value-of select="$sym"/></strong>        
                               </div>
                           </div>
                            <!-\- Now energies and dipole on two different tables -\->
                           <div class="row">
                               <div class="col-md-6 col-sm-12">                
                                   <xsl:variable name="energiesIndex" select="tokenize($energies/array[@dictRef='cc:serial'],'\s+')"></xsl:variable>
                                   <xsl:variable name="energy" select="tokenize($energies/array[@dictRef='cc:energy'],'\s+')"></xsl:variable>
                                   <xsl:variable name="oscillator" select="tokenize($energies/array[@dictRef='cc:oscillator'],'\s+')"></xsl:variable>
                                   <xsl:variable name="energyDiff" select="tokenize($energies/array[@dictRef='cc:energyDiff'],'\s+')"></xsl:variable>
                                   <span> Excitation energies E in a.u. , dE wrt prev. cycle,oscillator strengths f in a.u.</span>
                                   <table id="excitationEnergies-{generate-id($energies)}"></table>
                                   <script type="text/javascript">
                                       $(document).ready(function() {                        
                                            $("table#excitationEnergies-<xsl:value-of select="generate-id($energies)"/>").dataTable( { 
                                                "aaData" : [
                                                <xsl:for-each select="1 to count($energiesIndex)"> 
                                                    <xsl:variable name="innerIndex" select="."/>
                                                    [<xsl:value-of select="$energiesIndex[$innerIndex]"/>,'<xsl:value-of select="$energy[$innerIndex]"/>','<xsl:value-of select="$oscillator[$innerIndex]"/>','<xsl:value-of select="$energyDiff[$innerIndex]"/>']<xsl:if test="$innerIndex &lt; count($energiesIndex)"><xsl:text>,</xsl:text></xsl:if>                                                   
                                                </xsl:for-each>                                                                
                                                ],
                                                "aoColumns" : [
                                                    {"sTitle" : "no."},
                                                    {"sTitle" : "E/a.u."},
                                                    {"sTitle" : "f"},
                                                    {"sTitle" : "dE/a.u."}                                                    
                                                ],                                              
                                                "bFilter": false,                                 
                                                "bPaginate": false,                                    
                                                "bSort": false,
                                                "bInfo": false
                                             } );   
                                       } );                                       
                                   </script>
                               </div>
                               <div class="col-md-6 col-sm-12">
                                   <xsl:variable name="energiesIndex" select="tokenize($dipole/array[@dictRef='cc:serial'],'\s+')"></xsl:variable>
                                   <xsl:variable name="energy" select="tokenize($dipole/array[@dictRef='cc:energy'],'\s+')"></xsl:variable>
                                   <xsl:variable name="oscillator" select="tokenize($dipole/array[@dictRef='cc:oscillator' or @dictRef='cc:energyDiff'],'\s+')"></xsl:variable> <!-\- Bug fix -\->
                                   <xsl:variable name="muX" select="tokenize($dipole/array[@dictRef='cc:muX'],'\s+')"></xsl:variable>
                                   <xsl:variable name="muY" select="tokenize($dipole/array[@dictRef='cc:muY'],'\s+')"></xsl:variable>
                                   <xsl:variable name="muZ" select="tokenize($dipole/array[@dictRef='cc:muZ'],'\s+')"></xsl:variable>
                                   <span> Transition dipole moments mu (x,y,z) in a.u. (weak excitations are not printed)</span>                                   
                                   <table id="excitationDipole-{generate-id($dipole)}"></table>
                                   <script type="text/javascript">
                                       $(document).ready(function() {                        
                                       $("table#excitationDipole-<xsl:value-of select="generate-id($dipole)"/>").dataTable( { 
                                       "aaData" : [
                                       <xsl:for-each select="1 to count($energiesIndex)"> 
                                           <xsl:variable name="innerIndex" select="."/>
                                           [<xsl:value-of select="$energiesIndex[$innerIndex]"/>,
                                           '<xsl:value-of select="$energy[$innerIndex]"/>',
                                           '<xsl:value-of select="$oscillator[$innerIndex]"/>',
                                           '<xsl:value-of select="$muX[$innerIndex]"/>',
                                           '<xsl:value-of select="$muY[$innerIndex]"/>',
                                           '<xsl:value-of select="$muZ[$innerIndex]"/>']
                                           <xsl:if test="$innerIndex &lt; count($energiesIndex)"><xsl:text>,</xsl:text></xsl:if>                                                   
                                       </xsl:for-each>                                                                
                                       ],
                                       "aoColumns" : [
                                       {"sTitle" : "no."},
                                       {"sTitle" : "E/eV"},
                                       {"sTitle" : "f"},
                                       {"sTitle" : "X"},
                                       {"sTitle" : "Y"},
                                       {"sTitle" : "Z"}
                                       ],                                              
                                       "bFilter": false,                                 
                                       "bPaginate": false,                                    
                                       "bSort": false,
                                       "bInfo": false
                                       } );   
                                       } );                                       
                                   </script>
                               </div>
                           </div> 
                           

                           <xsl:variable name="excitationPlotId" select="generate-id($excitationEnergies[$outerIndex])" />
                            <div class="col-md-12">
                                <div id="excitationEnergyPlotlyContainer-{$excitationPlotId}" style="min-height:550px; width: 100%; display:none !important" class=" d-flex"/>
                            </div>                            
                            <xsl:variable name="energy" select="tokenize($energies/array[@dictRef='cc:energy'],'\s+')"></xsl:variable>
                            <xsl:variable name="oscillator" select="tokenize($energies/array[@dictRef='cc:oscillator'],'\s+')"></xsl:variable>
                            
                            <div class="col-sm-12 col-md-12">
                                <div id="excitationEnergyContainerControls-{$excitationPlotId}" class="mt-2 mb-5">                                
                                    <form>
                                        <div class="form-row align-items-center">
                                            <div class="col-auto">                                                
                                                <div class="form-group form-check-inline mt-4 mr-5">
                                                    <div id="excitationGraphUnits-{$excitationPlotId}" class="form-check">
                                                        <input class="form-check-input" type="radio" id="excitationGraphUnits-{$excitationPlotId}-ev" name="excitationGraphUnits-{$excitationPlotId}" value="ev" onclick="excitationParametersChanged('{$excitationPlotId}')" />
                                                        <label class="form-check-label" for="excitationGraphUnits-{$excitationPlotId}-ev">eV</label>                                                                        
                                                    </div>    
                                                    <div class="form-check">
                                                        <input class="form-check-input" type="radio" id="excitationGraphUnits-{$excitationPlotId}-cm1" name="excitationGraphUnits-{$excitationPlotId}" value="cm1" onclick="excitationParametersChanged('{$excitationPlotId}')" />
                                                        <label class="form-check-label" for="excitationGraphUnits-{$excitationPlotId}-cm1">cm-1</label>                                                                        
                                                    </div>
                                                    <div class="form-check">
                                                        <input class="form-check-input" type="radio" id="excitationGraphUnits-{$excitationPlotId}-nm" name="excitationGraphUnits-{$excitationPlotId}" value="nm" onclick="excitationParametersChanged('{$excitationPlotId}')" />
                                                        <label class="form-check-label" for="excitationGraphUnits-{$excitationPlotId}-nm">nm</label>                                                                        
                                                    </div>
                                                </div>
                                            </div>
                                            <div class="col-auto">
                                                <div class="form-group">
                                                    <label for="bandwidth-{$excitationPlotId}">Bandwidth</label>
                                                    <input type="text" class="form-control" id="bandwidth-{$excitationPlotId}" value="0.15"/>
                                                </div>
                                            </div>
                                            <div class="col-auto">
                                                <div class="form-group">
                                                    <label for="minXEnergy-{$excitationPlotId}">min X</label>
                                                    <input type="text" class="form-control" id="minXEnergy-{$excitationPlotId}" value="0"/>
                                                </div>
                                            </div>
                                            <div class="col-auto">
                                                <div class="form-group">
                                                    <label for="maxXEnergy-{$excitationPlotId}">max X</label>
                                                    <input type="text" class="form-control" id="maxXEnergy-{$excitationPlotId}" value="10"/>
                                                </div>
                                            </div>
                                            <div class="col-auto">
                                                <div class="form-group mt-4">
                                                    <button type="button" class="btn btn-secondary" disabled="disabled" id="excitation-generateBtn-{$excitationPlotId}" onclick="generateExcitationChart('{$excitationPlotId}')">Plot!</button>
                                                </div>
                                            </div>
                                        </div>
                                    </form>
                                </div>
                            </div>                           
                            
                            <xsl:variable name="oscStrengths" as="xs:float*">                               
                                <xsl:for-each select="$oscillator"><xsl:element name="value"><xsl:value-of select="."/></xsl:element></xsl:for-each>
                            </xsl:variable>
                            <xsl:variable name="energyEv" as="xs:float*">                               
                                <xsl:for-each select="$energy"><xsl:element name="value"><xsl:value-of select="(number(.)*$AU_TO_EV)"/></xsl:element></xsl:for-each>                                                   
                            </xsl:variable>
                            <xsl:variable name="energyCm1" as="xs:float*">                             
                                <xsl:for-each select="$energy"><xsl:element name="value"><xsl:value-of select="(number(.)*$AU_TO_EV) * $EV_TO_CM_1"/></xsl:element></xsl:for-each>
                            </xsl:variable>
                            <xsl:variable name="energyNm" as="xs:float*">                              
                                <xsl:for-each select="$energy"><xsl:element name="value"><xsl:value-of select="format-number($EV_TO_NM div (number(.) * $AU_TO_EV),'#0.000')"/></xsl:element></xsl:for-each>
                            </xsl:variable>
                            <script type="text/javascript">                                                                  
                                oscStrengths['<xsl:value-of select="$excitationPlotId"/>'] = [<xsl:value-of select="helper:toTextArray($oscStrengths)"/>];                        
                                energyEv['<xsl:value-of select="$excitationPlotId"/>'] = [<xsl:value-of select="helper:toTextArray($energyEv)"/>];                     
                                energyCm1['<xsl:value-of select="$excitationPlotId"/>'] = [<xsl:value-of select="helper:toTextArray($energyCm1)"/>];
                                energyNm['<xsl:value-of select="$excitationPlotId"/>'] = [<xsl:value-of select="helper:toTextArray($energyNm)"/>];
                                
                                <xsl:variable name="energyMax" select="number(format-number( round(10*max($energyEv)) div 10 ,'##0.0'))"/>
                                <xsl:variable name="energyMin" select="number(format-number( round(10*min($energyEv)) div 10 ,'##0.0'))"/>
                                
                                energyMaxEv['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$energyMax"/> + evBandwidth);                 
                                energyMinEv['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$energyMin"/> - evBandwidth);                
                                energyMaxCm1['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$energyMax * $EV_TO_CM_1"/> + cm1Bandwidth);          
                                energyMinCm1['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$energyMin * $EV_TO_CM_1"/> - cm1Bandwidth);
                                energyMaxNm['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$EV_TO_NM div $energyMin"/> + nmBandwidth);   
                                energyMinNm['<xsl:value-of select="$excitationPlotId"/>'] = round(<xsl:value-of select="$EV_TO_NM div $energyMax"/> - nmBandwidth);
                            </script>                                                
                        </xsl:for-each>
                    </div>
                </div>
            </div>
        </xsl:if>        
    </xsl:template>-->
    
 <!--   <!-\- Rotatory strengths -\->
    <xsl:template name="rotatorystrengths">
        <xsl:variable name="excitations" select=".//module[@id='finalization']//module[@cmlx:templateRef='excitation.energy']"/>
        <xsl:if test="exists($excitations//module[@cmlx:templateRef='rotatory'])">
            <xsl:variable name="excitationEnergies" select="$excitations/module[@cmlx:templateRef='excitationEnergies']"/>             
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#excitations-{generate-id($excitations)}-rotatory" style="cursor: pointer; cursor: hand;">                    
                    <h4 class="panel-title">
                        Rotatory strengths
                    </h4>
                </div>
                <div id="excitations-{generate-id($excitations)}-rotatory" class="panel-collapse collapse">
                    <div class="panel-body">    
                        <xsl:for-each select="1 to count($excitationEnergies)">
                            <xsl:variable name="outerIndex" select="."/>
                            <xsl:variable name="energies"       select="$excitationEnergies[$outerIndex]/module[@cmlx:templateRef='energies']"/>                            
                            <xsl:variable name="rotstrengths"   select="$excitationEnergies[$outerIndex]/module[@cmlx:templateRef='rotatory']"/>
                            <xsl:variable name="sym" select="$excitationEnergies[$outerIndex]/scalar[@dictRef='cc:symm']"></xsl:variable>
                            <!-\- Print header -\->
                            <div class="row">
                                <div class="col-md-12">
                                    <strong>Symmetry <xsl:value-of select="$sym"/></strong>        
                                </div>
                            </div>
                            <!-\- Now energies and rotatory strengths on two different tables -\->
                            <div class="row">
                                <div class="col-md-6 col-sm-12">                
                                    <xsl:variable name="energiesIndex" select="tokenize($energies/array[@dictRef='cc:serial'],'\s+')"></xsl:variable>
                                    <xsl:variable name="energy" select="tokenize($energies/array[@dictRef='cc:energy'],'\s+')"></xsl:variable>
                                    <xsl:variable name="oscillator" select="tokenize($energies/array[@dictRef='cc:oscillator'],'\s+')"></xsl:variable>
                                    <xsl:variable name="energyDiff" select="tokenize($energies/array[@dictRef='cc:energyDiff'],'\s+')"></xsl:variable>
                                    <span> Excitation energies E in a.u. , dE wrt prev. cycle,oscillator strengths f in a.u.</span><br/><br/>
                                    <table id="excitationEnergies-{generate-id($energies)}-rotatory"></table>
                                    <script type="text/javascript">
                                        $(document).ready(function() {                        
                                        $("table#excitationEnergies-<xsl:value-of select="generate-id($energies)"/>-rotatory").dataTable( { 
                                        "aaData" : [
                                        <xsl:for-each select="1 to count($energiesIndex)"> 
                                            <xsl:variable name="innerIndex" select="."/>
                                            [<xsl:value-of select="$energiesIndex[$innerIndex]"/>,'<xsl:value-of select="$energy[$innerIndex]"/>','<xsl:value-of select="$oscillator[$innerIndex]"/>','<xsl:value-of select="$energyDiff[$innerIndex]"/>']<xsl:if test="$innerIndex &lt; count($energiesIndex)"><xsl:text>,</xsl:text></xsl:if>                                                   
                                        </xsl:for-each>                                                                
                                        ],
                                        "aoColumns" : [
                                        {"sTitle" : "no."},
                                        {"sTitle" : "E/a.u."},
                                        {"sTitle" : "f"},
                                        {"sTitle" : "dE/a.u."}                                                    
                                        ],                                              
                                        "bFilter": false,                                 
                                        "bPaginate": false,                                    
                                        "bSort": false,
                                        "bInfo": false
                                        } );   
                                        } );                                       
                                    </script>
                                </div>
                                <div class="col-md-6 col-sm-12">
                                    <xsl:variable name="rotatoryIndex" select="tokenize($rotstrengths/array[@dictRef='cc:serial'],'\s+')"></xsl:variable>
                                    <xsl:variable name="strengths" select="tokenize($rotstrengths/array[@dictRef='a:strengths'],'\s+')"></xsl:variable>
                                    <xsl:variable name="mX" select="tokenize($rotstrengths/array[@dictRef='cc:mX'],'\s+')"></xsl:variable>
                                    <xsl:variable name="mY" select="tokenize($rotstrengths/array[@dictRef='cc:mY'],'\s+')"></xsl:variable>
                                    <xsl:variable name="mZ" select="tokenize($rotstrengths/array[@dictRef='cc:mZ'],'\s+')"></xsl:variable>
                                    <span>  Rotatory strengths R in 10**(-40) esu**2 * cm**2 , (multiply by 1.07827 to obtain reduced rotatory strengths), magnetic transition dipole vectors m in a.u.:</span>                                   
                                    <table id="excitationRotationStregths-{generate-id($rotstrengths)}"></table>
                                    <script type="text/javascript">
                                        $(document).ready(function() {                        
                                        $("table#excitationRotationStregths-<xsl:value-of select="generate-id($rotstrengths)"/>").dataTable( { 
                                        "aaData" : [
                                        <xsl:for-each select="1 to count($rotatoryIndex)"> 
                                            <xsl:variable name="innerIndex" select="."/>
                                            [<xsl:value-of select="$rotatoryIndex[$innerIndex]"/>,
                                            '<xsl:value-of select="$strengths[$innerIndex]"/>',                                            
                                            '<xsl:value-of select="$mX[$innerIndex]"/>',
                                            '<xsl:value-of select="$mY[$innerIndex]"/>',
                                            '<xsl:value-of select="$mZ[$innerIndex]"/>']
                                            <xsl:if test="$innerIndex &lt; count($rotatoryIndex)"><xsl:text>,</xsl:text></xsl:if>                                                   
                                        </xsl:for-each>                                                                
                                        ],
                                        "aoColumns" : [
                                        {"sTitle" : "no."},
                                        {"sTitle" : "R"},
                                        {"sTitle" : "X"},
                                        {"sTitle" : "Y"},
                                        {"sTitle" : "Z"}
                                        ],                                              
                                        "bFilter": false,                                 
                                        "bPaginate": false,                                    
                                        "bSort": false,
                                        "bInfo": false
                                        } );   
                                        } );                                       
                                    </script>
                                </div>
                            </div>                                                       
                        </xsl:for-each>
                    </div>
                </div>
            </div>
        </xsl:if>        
    </xsl:template>-->
    
    <xsl:function name="helper:toTextArray">
        <xsl:param name="sequence" as="xs:float*"/>
        <xsl:for-each select="1 to count($sequence)">
            <xsl:variable name="outerIndex" select="."/>
            <xsl:value-of select="$sequence[$outerIndex]"/>
            <xsl:if test="$outerIndex &lt; count($sequence)">, </xsl:if>
        </xsl:for-each>        
    </xsl:function>
 <!--   <!-\- Print license footer -\->
    <xsl:template name="printLicense">
        <div class="row">
            <div class="col-md-12 text-right">
                <br/>
                <span>Report data <a rel="license" href="http://creativecommons.org/licenses/by/4.0/" target="_blank"><img alt="Creative Commons License" style="border-width:0" src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAMAAABEF7i9AAAABGdBTUEAANbY1E9YMgAAAJZQTFRF////7u7u3d3dys7KzMzMyMzIxsrGxcbFur+6u7u7s7iyq7GqqqqqmZmZlJmTj5CPiIiIh4eHhoaGgICAfYJ9d3d3cnZxZ2tnZmZmW15bVVVVS0xLREREQ0NDQkJCQUJBOz07OTs5MzMzMTMxLjAuJygnJCUjIiIiISEhICAgGRkZERERDxAPDg4ODQ4NDQ0NDQ0MAAAADbeuvgAAAMNJREFUeNqtk9sSgiAQQJekslaNLpgZ3exiZWr7/z/XEI6N00spO/DCMocDywJZDiCwGhqIOpYkqiWVB4jOo7WfAQa5qA9RZ0R33hG4v36sOa0Q++tuwEIAT0kxRSkHdUR0Ju88gN5k5j/ABY0gjUHKjPkeyALRHZq8GfCvYUic6arEib6zzBHHg9rwd78vQxVlTPogy6ZhCyCWAryMIqYo6cHm1HjDVsDDzXKVg+e0Bm4vFv4hhjSreLt7106x3cuW4wXCCZ5As6hO5AAAAABJRU5ErkJggg==" /></a></span>   
                <br/>
                <span>This HTML file <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/" target="_blank"><img alt="Creative Commons License" style="border-width:0" src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAMAAABEF7i9AAAABGdBTUEAANbY1E9YMgAAAJZQTFRF////7u7u3d3dys7KzMzMyMzIxsrGxcbFur+6u7u7s7iyq7GqqqqqmZmZlJmTj5CPiIiIh4eHhoaGgICAfYJ9d3d3cnZxZ2tnZmZmW15bVVVVS0xLREREQ0NDQkJCQUJBOz07OTs5MzMzMTMxLjAuJygnJCUjIiIiISEhICAgGRkZERERDxAPDg4ODQ4NDQ0NDQ0MAAAADbeuvgAAAOJJREFUeNqtk21PAjEQhKdyKjpwFrGIWN8QOPE8cP7/n/PD9YolmhivmzRpJ9un090WyhwQsoYgkCRvZYPkm5IcfPzbXwssGxsP8WtyOO0JfH47uC50R57e9wPuLIpK3nhVBfwrObiSJAAS1A7FKXAAhDYcAW90gWoB52ozHsHtyOF5yEebFfK71W9CB5ypMLLAYgkAriEvz6LD4KPb2+FTAT869PPauDHcPnWo7zdE6vBIiDXcW4xqzY3X8gQPq6SGKfBvNeTLNnOXy89JBD5uMrxDznQdeE/vfX9K7r+cOb4AY2+UGwcd6o0AAAAASUVORK5CYII=" /></a></span>                
            </div>            
        </div>
    </xsl:template>    -->
        
    <!-- Override default templates -->
    <xsl:template match="text()"/>
    <xsl:template match="*"/>    
	
</xsl:stylesheet>
