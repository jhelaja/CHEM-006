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
    xmlns:vasp="https://www.vasp.at/"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"

    xpath-default-namespace="http://www.xml-cml.org/schema"
    exclude-result-prefixes="xs xd cml ckbk vasp helper cmlx"
    version="2.0">

    <xd:doc scope="stylesheet">
        <xd:desc>
            <xd:p><xd:b>Created on:</xd:b>May 30 2024</xd:p>
            <xd:p><xd:b>Author:</xd:b>Moisés Álvarez Moreno and Diego Garay Ruiz</xd:p>
            <xd:p><xd:b>Center:</xd:b>Institute of Chemical Research of Catalonia</xd:p>            
        </xd:desc>
    </xd:doc>
    <xsl:include href="vasp_helper.xsl"/>
    <xsl:include href="chem_helper.xsl"/>
    <xsl:output method="text" omit-xml-declaration="yes" indent="no" />
    <xsl:strip-space elements="*"/>


    <xsl:param name="title"/>
    <xsl:param name="author"/>

    <xsl:variable name="ibrion" select="//module[@id='job']/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ibrion']/scalar"/>    
    <xsl:variable name="calcType" select="vasp:getCalcType($ibrion)"/>    
    <xsl:variable name="ispin" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ispin']/scalar"/>
    <xsl:variable name="shellType" select="vasp:getShellType($ispin)"/>
    <xsl:variable name="ediff" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ediff']/scalar"/>
    <xsl:variable name="ediffg" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ediffg']/scalar"/>
    <xsl:variable name="potim" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:potim']/scalar"/>
    <xsl:variable name="ldipol" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ldipol']/scalar"/>
    <xsl:variable name="idipol" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:idipol']/scalar"/>
    <xsl:variable name="grimmes" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/module[@id='otherComponents']/module[@cmlx:templateRef='grimmes']"/>

    <xsl:variable name="energyCutoff" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:energyCutoff']/scalar"/>
    
    <xsl:variable name="isNeb" select="count($ibrion) > 1 and not(exists($ibrion[text() = '44']))"/>

    <!-- Functional calculation -->   
    <xsl:variable name="gga"      select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:gga']/scalar"/>
    <xsl:variable name="lhfcalc"  select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:lhfcalc']/scalar"/>
    <xsl:variable name="ldau"     select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ldau']/scalar"/>
    <xsl:variable name="aggac"    select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:aggac']/scalar"/>
    <xsl:variable name="hfscreen" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:hfscreen']/scalar"/>    
    <xsl:variable name="luseVdw"  select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:luseVdw']/scalar"/>
    <xsl:variable name="zabVdw"   select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:zabVdw']/scalar"/>
    <xsl:variable name="param1"   select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:param1']/scalar"/>    
    <xsl:variable name="param2"   select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:param2']/scalar"/>    
    <xsl:variable name="aexx" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:aexx']/scalar"/>
    <xsl:variable name="aggax" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:aggax']/scalar"/>
    <xsl:variable name="aldac" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:aldac']/scalar"/>
    <xsl:variable name="functional" select="vasp:getFunctional($gga,$lhfcalc,$aggac,$hfscreen,$luseVdw,$zabVdw, $param1, $param2, $ldau, $aexx, $aggax, $aldac)"/>
    
    <!-- Environment module -->
    <xsl:variable name="programParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:program']"/>
    <xsl:variable name="versionParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:programVersion']"/>
    <!-- Initialization module -->    
    <xsl:variable name="multiplicity" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:nupdown']/scalar"/>
    <xsl:variable name="sigma"        select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:sigma']/scalar"/>
    <xsl:variable name="ismear"       select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ismear']/scalar"/>
    <xsl:variable name="nelect"       select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:nelect']/scalar"/>
    <xsl:variable name="kpoints"      select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='kpoints']"/>
    <xsl:variable name="potcar"       select="(//module[@id='job'][1]/module[@dictRef='cc:initialization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='potcar'])[1]"/>
    <xsl:variable name="temperature" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='cc:temp']/scalar"/>
    <!-- Calculation module -->
    <xsl:variable name="forces" select="//module[@id='job'][1]/module[@dictRef='cc:calculation']/module[@dictRef='cc:userDefinedModule']/module[@id='forces']"/>
    <xsl:variable name="dos" select="//module[@id='job'][1]/module[@dictRef='cc:calculation']/module[@dictRef='cc:userDefinedModule']/module[@id='vasp.doscar']"/>    
    <!-- Finalization module -->
    <xsl:variable name="magnetizations" select="//module[@dictRef='cc:finalization']/propertyList/property/module[@cmlx:templateRef='magnetization']"/>
    <xsl:variable name="vibrations" select="(//module[@dictRef='cc:finalization']/propertyList/property/module[@cmlx:templateRef='vibrations'])[last()]"/>        
    <!-- Geometry -->
    <xsl:variable name="energies" select="//module[@id='job']/module[@dictRef='cc:finalization']//property[@dictRef='cc:energies']/module[@cmlx:templateRef='energies']"/>
    <xsl:variable name="initialMolecule" select="//module[@dictRef='cc:initialization']/molecule[@id='initial']" />
    <xsl:variable name="isOptimization" select="compare($calcType, $vasp:GeometryOptimization) = 0" />
    <xsl:variable name="finalMolecules" select="//module[@dictRef='cc:finalization' and child::molecule]//molecule"/>    
    <xsl:variable name="isLargeMolecule" select="count($finalMolecules[1]/atomArray/atom) &gt; 1000" />

    <xsl:template match="/">
        <xsl:call-template name="generalInfo"/>
        <xsl:call-template name="settings"/>
        <xsl:choose>
            <xsl:when test="count($finalMolecules) = 1">
                <xsl:call-template name="atomicCoordinatesSingle">
                    <xsl:with-param name="molecule" select="$finalMolecules"/>                                                                                                                          
                </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
                <xsl:call-template name="atomicCoordinatesMultiple">
                    <xsl:with-param name="molecules" select="$finalMolecules"/>                                                                                                                          
                </xsl:call-template> 
            </xsl:otherwise>
        </xsl:choose>
        <!-- Molecular Info -->                    
        <xsl:if test="exists($kpoints)">                                                  
            <xsl:call-template name="kpoints"/>
        </xsl:if>                                        
        <xsl:call-template name="energy"/>
<!--        <xsl:call-template name="eigenvalues"/>-->

        <xsl:call-template name="magnetization"/>
        <xsl:call-template name="vibrations"/>
        <xsl-text>GENERAL_END&#xa;</xsl-text>

    </xsl:template>
    
    <!-- General Info -->
    <xsl:template name="generalInfo">
            <xsl:if test="$title">
                <xsl:value-of select="concat('title:',$title,'#;#')"/>
            </xsl:if>
            <xsl:value-of select="concat('program:',$programParameter/scalar/text(),'#;#')"/>                                        
            <xsl:value-of select="concat('version:',$versionParameter/scalar/text(),'#;#')"/>
            <xsl:if test="$author">
                <xsl:value-of select="concat('author:',$author,'#;#')"/>
            </xsl:if>
            <xsl:variable name="calcFormula" select="helper:calculateHillNotationFormula($finalMolecules[1]/atomArray/atom)"/>
            <xsl:value-of select="concat('formula:',string-join($calcFormula),'#;#')"/>

            <xsl:value-of select="concat('calcType:',$calcType,'#;#')"/>                                        
            <xsl:value-of select="concat('functional:',($functional)[1],'#;#')"/>
            <xsl:value-of select="concat('shellType:',$shellType,'#;#','ispin:',$ispin,'#;#')"/>                    
            
            <xsl:if test="exists($temperature)">                                  
                <xsl:value-of select="concat('temperature:',$temperature/text(),'#;#')"></xsl:value-of>   
                <xsl:value-of select="concat('temperatureUnits:',helper:printUnitSymbol($temperature/@units),'#;#')"/>
            </xsl:if>            
    </xsl:template>
    
    <xsl:template name="settings">
        <xsl:value-of select="concat('sigma:',$sigma,'#;#')"/>
        <xsl:value-of select="concat('ismear:',$ismear,'#;#')"/>
        <xsl:if test="exists($ldipol)">
            <xsl:value-of select="concat('ldipol:',$ldipol,'#;#')"/>                                         
        </xsl:if>
        <xsl:if test="exists($idipol)">
            <xsl:value-of select="concat('idipol:',$idipol,'#;#')"/>
        </xsl:if>
        
        <xsl:if test="exists($multiplicity) and $multiplicity &gt; -1">
            <xsl:value-of select="concat('multiplicity:',$multiplicity,'#;#')"/>                                
        </xsl:if>
        <xsl:value-of select="concat('nelect:',$nelect,'#;#')"/>                            
        <xsl:if test="exists($energyCutoff)">
            <xsl:value-of select="concat('energyCutoff:',format-number($energyCutoff,'#0.00'),'#;#')"/>
        </xsl:if>
        <xsl:value-of select="concat('ediff:',$ediff,'#;#')"/>
        <xsl:if test="not(($calcType = $vasp:FrequencyCalculus ) or ( $calcType = $vasp:SinglePoint)) ">
            <xsl:value-of select="concat('ediffg:',$ediffg,'#;#')"/>
        </xsl:if>
        <xsl:value-of select="concat('potim:',$potim,'#;#')"/>
        <xsl:if test="exists(//parameter[contains(@dictRef,'v:ldau')])">
            <xsl:variable name="ldaul" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ldaul']/array"/>
            <xsl:variable name="ldauu" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ldauu']/array"/>
            <xsl:variable name="ldauj" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='v:ldauj']/array"/>
            <xsl:value-of select="concat('ldaul:',$ldaul,'#;#')" />
            <xsl:value-of select="concat('ldauu:',$ldauu,'#;#')" />
            <xsl:value-of select="concat('ldauj:',$ldauj,'#;#')" />
        </xsl:if>
        <xsl:if test="exists(//parameter[contains(@dictRef,'v:lvdw')]) and compare(//parameter[contains(@dictRef,'v:lvdw')]/scalar,'true') = 0">
            <xsl-text>lvdw:T#;#</xsl-text>
            <xsl:choose>
                <xsl:when test="not(exists(//parameter[contains(@dictRef,'v:vdwversion')]))">
                    <xsl-text>vdwVersion:D2#;#</xsl-text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl:value-of select="concat('vdwVersion:D',//parameter[contains(@dictRef,'v:vdwversion')],'#;#')"/>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:if>
        
        <xsl:if test="exists($grimmes)">  
            <xsl:variable name="atomType" select="tokenize($grimmes/array[@dictRef='cc:elementType'],'\s+')"/>
            <xsl:variable name="grimmeC6" select="tokenize($grimmes/array[@dictRef='v:grimmeC6'],'\s+')"/>
            <xsl:variable name="grimmeR0" select="tokenize($grimmes/array[@dictRef='v:grimmeR0'],'\s+')"/>
            <xsl-text>grimmeParameters:atom C6 R0 '&#xa;'</xsl-text>
            <xsl:for-each select="1 to count($atomType)">
                <xsl:variable name="outerIndex" select="."/>
                <xsl:value-of select="$atomType[$outerIndex]"/>
                <xsl:value-of select="$grimmeC6[$outerIndex]"/>
                <xsl:value-of select="$grimmeR0[$outerIndex]"/>
                <xsl-text>'&#xa;'</xsl-text>
            </xsl:for-each>      
            <xsl-text>grimmeParameterUnits:None Jnm^6/mol A#;#</xsl-text>   
        </xsl:if>
        
    </xsl:template>
    
    
    <!-- Atomic coordinates -->
    <xsl:template name="atomicCoordinatesSingle">
        <xsl:param name="molecule"/>

        <xsl:variable name="basisPerElement">
            <xsl:variable name="elements" select="tokenize($potcar/array[@dictRef = 'cc:atomType'], '\s+')"/>
            <xsl:variable name="atomsPerType" select="tokenize($potcar/array[@dictRef = 'cc:atomcount'], '\s+')"/>
            <xsl:variable name="basis" select="tokenize(replace($potcar/array[@dictRef = 'v:pseudopotential'], '\n', ''), '\|')"/>                                                                               
            <xsl:for-each select="1 to count($elements)">
                <xsl:variable name="outerIndex" select="position()"/>
            </xsl:for-each>      
        </xsl:variable>                                                                

        <!-- multiple loops to show cartesian, fractional coordinates, and basis-->
        <xsl-text>geometryCartesian:</xsl-text>
        <xsl:for-each select="$molecule/cml:atomArray/cml:atom">
            <xsl:variable name="elementType" select="./@elementType"/>                                                      
            <xsl:value-of select="concat($elementType,' ',format-number(@x3, '#0.0000'),' ',format-number(@y3, '#0.0000'),' ',format-number(@z3, '#0.0000'),'&#xa;')"/>
        </xsl:for-each> 
        <xsl-text>#;#</xsl-text>
        <xsl-text>geometryFractional:</xsl-text>
        <xsl:for-each select="$molecule/cml:atomArray/cml:atom">
            <xsl:variable name="elementType" select="./@elementType"/>                                                      
            <xsl:value-of select="concat($elementType,' ',format-number(@xFract, '#0.0000'),' ',format-number(@yFract, '#0.0000'),' ',format-number(@zFract, '#0.0000'),'&#xa;')"/>
        </xsl:for-each> 
        <xsl-text>#;#</xsl-text>
        <xsl-text> speciesBasis: </xsl-text>
        <xsl:for-each select="$molecule/cml:atomArray/cml:atom">
            <xsl:variable name="elementType" select="./@elementType"/>                                                      
            <xsl:variable name="id" select="@id"/>
            <xsl:variable name="outerIndex" select="position()"/>
            <xsl:variable name="basis" select="$basisPerElement/node()/node()[$outerIndex]/text()"/>
            <xsl:value-of select="concat($elementType,' ',$basis)"/>
        </xsl:for-each>
        <xsl-text>#;#</xsl-text>
        <xsl-text>cellParameters:</xsl-text>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='a'],' ')"/>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='b'],' ')"/>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='c'],' ')"/>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='alpha'],' ')"/>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='beta'],' ')"/>
        <xsl:value-of select="concat($molecule/cml:crystal/scalar[@title='gamma'],' ')"/>
        <xsl-text>#;#</xsl-text>
    </xsl:template>
  
    <xsl:template name="atomicCoordinatesMultiple">
        <xsl:param name="molecules"/>

<!--        <xsl:apply-templates select="$finalMolecules[1]/crystal" />
        <xsl:apply-templates select="$potcar/array[@dictRef='cc:valence']" />-->

<!--        <xsl:value-of select="count(./cml:atomArray/cml:atom)"/>-->
        <xsl:variable name="basisPerElement">
                <xsl:variable name="elements" select="tokenize($potcar/array[@dictRef = 'cc:atomType'], '\s+')"/>
                <xsl:variable name="atomsPerType" select="tokenize($potcar/array[@dictRef = 'cc:atomcount'], '\s+')"/>
                <xsl:variable name="basis" select="tokenize(replace($potcar/array[@dictRef = 'v:pseudopotential'], '\n', ''), '\|')"/>                                                                               
                <xsl:for-each select="1 to count($elements)">
                    <xsl:variable name="outerIndex" select="position()"/>
                    <!--<xsl:for-each select="1 to xs:integer($atomsPerType[$outerIndex])">                                                                                        
                        <element type="{$elements[$outerIndex]}"><xsl:value-of select="$basis[$outerIndex]"/></element>
                    </xsl:for-each>-->
                </xsl:for-each>      
        </xsl:variable>                                                                
        <xsl:for-each select="$molecules">
            <xsl:variable name="molecule" select="."/>
            <xsl:variable name="index" select="position()"/>
            
            <!-- multiple loops to show cartesian, fractional coordinates, and basis-->

            <xsl:value-of select="concat('geometryCartesian',$index,':')"/>
            <xsl:for-each select="./cml:atomArray/cml:atom">
                    <xsl:variable name="elementType" select="./@elementType"/>                                                      
                    <xsl:value-of select="concat($elementType,' ',format-number(@x3, '#0.0000'),' ',format-number(@y3, '#0.0000'),' ',format-number(@z3, '#0.0000'),'&#xa;')"/>
            </xsl:for-each> 
            <xsl-text>#;#</xsl-text>
            <xsl:value-of select="concat('geometryFractional',$index,':')"/>
            <xsl:for-each select="./cml:atomArray/cml:atom">
                <xsl:variable name="elementType" select="./@elementType"/>                                                      
                <xsl:value-of select="concat($elementType,' ',format-number(@xFract, '#0.0000'),' ',format-number(@yFract, '#0.0000'),' ',format-number(@zFract, '#0.0000'),'&#xa;')"/>
            </xsl:for-each> 
            <xsl-text>#;#</xsl-text>
            <xsl:value-of select="concat('speciesBasis',$index,':')"/>
            <xsl:for-each select="./cml:atomArray/cml:atom">
                <xsl:variable name="elementType" select="./@elementType"/>                                                      
                <xsl:variable name="id" select="@id"/>
                <xsl:variable name="outerIndex" select="position()"/>
                <xsl:variable name="basis" select="$basisPerElement/node()/node()[$outerIndex]/text()"/>
                <xsl:value-of select="concat($elementType,' ',$basis)"/>
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>
            <xsl:value-of select="concat('cellParameters',$index,':')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='a'],' ')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='b'],' ')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='c'],' ')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='alpha'],' ')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='beta'],' ')"/>
            <xsl:value-of select="concat(./cml:crystal/scalar[@title='gamma'],' ')"/>
            <xsl-text>#;#</xsl-text>

        </xsl:for-each>
    </xsl:template>

    <!-- kpoint list -->
    <xsl:template name="kpoints">
        <xsl:variable name="subdivisionN" select="tokenize($kpoints/array[@dictRef='v:subdivisionN'],'\s+')"/>
        <xsl:variable name="shiftS" select="tokenize($kpoints/array[@dictRef='v:shiftS'],'\s+')"/>
        
                <xsl:if test="exists($subdivisionN)">
                    <xsl-text>kpointSubdivisions: </xsl-text>
                     <xsl:for-each select="1 to count($subdivisionN)">
                         <xsl:variable name="outerIndex" select="."/>
                         <xsl:value-of select="$subdivisionN[$outerIndex]"/>
                     </xsl:for-each>
                    <xsl-text>#;#</xsl-text>
                    <xsl-text>kpointShifts: </xsl-text>
                    <xsl:for-each select="1 to count($shiftS)">
                        <xsl:variable name="outerIndex" select="."/>
                        <xsl:value-of select="$shiftS[$outerIndex]"/>
                    </xsl:for-each>
                    <xsl-text>#;#</xsl-text>
                </xsl:if>

            <xsl:if test="exists($kpoints/array[@dictRef='v:kpointlist'])">
                <xsl:variable name="coords" select="tokenize($kpoints/array[@dictRef='v:kpointlist'], '\s+')"></xsl:variable>                    
                <xsl:variable name="weights" select="tokenize($kpoints/array[@dictRef='v:weight'], '\s+')"></xsl:variable>
                <xsl-text>coordinateSpace:</xsl-text>                  
                        <xsl:choose>
                            <xsl:when test="$kpoints/scalar[@dictRef='v:coordtype' and matches(lower-case(text()), 'cartesian')]">Cartesian </xsl:when>
                            <xsl:otherwise>Reciprocal </xsl:otherwise>
                        </xsl:choose>
                <xsl-text>#;#</xsl-text>  
                <xsl:for-each select="1 to xs:integer(number($kpoints/array[@dictRef='v:kpointlist']/@size) div 3)">
                    <xsl:variable name="outerIndex" select="."/>
                    <xsl:variable name="kpointx" select="$coords[3 * ($outerIndex - 1) + 1]"/>
                    <xsl:variable name="kpointy" select="$coords[3 * ($outerIndex - 1) + 2]"/>
                    <xsl:variable name="kpointz" select="$coords[3 * ($outerIndex - 1) + 3]"/>
                    <xsl:value-of select="concat($kpointx,' ',$kpointy,' ',$kpointz,' ')"/>
                    <xsl:choose>
                        <xsl:when test="exists($weights)">
                            <xsl:value-of select="concat($weights[$outerIndex],'&#xa;')"/>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl-text> &#xa; </xsl-text>/>
                        </xsl:otherwise>
                    </xsl:choose>
                </xsl:for-each>
            </xsl:if>                  
    </xsl:template>    
    
    <!-- Final energy -->
    <xsl:template name="energy">        
        <xsl:if test="count($energies) = 1">    <!-- Single energy field, display all energies -->
            <xsl:choose>
                <xsl:when test="$calcType != $vasp:MolecularDynamics">
                    <xsl:variable name="energy" select="$energies[1]"/>
                    <xsl:value-of select="concat('gibbsFreeEnergy:',$energy/*[@dictRef='cc:freeEnergy'],'#;#')"/>
                    <xsl:value-of select="concat('gibbsFreeEnergyUnits:',helper:printUnitSymbol($energy/*[@dictRef='cc:freeEnergy']/@units),'#;#')"/>
                    <xsl:value-of select="concat('e0:',$energy/*[@dictRef='cc:e0'],'#;#')"/>
                    <xsl:value-of select="concat('e0Units:',helper:printUnitSymbol($energy/*[@dictRef='cc:e0']/@units),'#;#')"/>

                    <xsl:if test="exists($energy/*[@dictRef='cc:deltaEnergy'])">
                            <xsl:value-of select="concat('deltaE:',number($energy/*[@dictRef='cc:deltaEnergy']),'#;#')"/>
                            <xsl:value-of select="concat('deltaEUnits:',helper:printUnitSymbol($energy/*[@dictRef='cc:deltaEnergy']/@units),'#;#')"/>
                    </xsl:if>                        
                    <xsl:if test="exists($energy/*[@dictRef='v:efermi'])">
                        <xsl:value-of select="concat('Efermi:',number($energy/*[@dictRef='v:efermi']),'#;#')"/>                    
                        <xsl:value-of select="concat('EfermiUnits:',helper:printUnitSymbol($energy/*[@dictRef='v:efermi']/@units),'#;#')"/>
                    </xsl:if>        
                </xsl:when>
                <xsl:otherwise>                                                
                    <xsl:variable name="free" select="tokenize($energies/cml:array[@dictRef='cc:freeEnergy'], '\s+')" />
                    <xsl:variable name="noentropy" select="tokenize($energies/cml:array[@dictRef='v:noEntropyEnergy'], '\s+')" />
                    <xsl:variable name="e0" select="tokenize($energies/cml:array[@dictRef='cc:e0'], '\s+')" />
                    <xsl:variable name="fermi" select="tokenize($energies/cml:array[@dictRef='v:efermi'], '\s+')" />
                    <xsl:variable name="delta" select="tokenize(concat('0.0 ', $energies/cml:array[@dictRef='cc:deltaEnergy']), '\s+')" />
                    <xsl:variable name="f" select="tokenize($forces/cml:matrix[@dictRef='cc:force'], '\s+')" />
                    <xsl:variable name="atomCount" select="xs:integer(number($forces/cml:matrix[@dictRef='cc:force']/@cols) div 3)" />
                    <xsl:variable name="maxforces">
                        <xsl:for-each select="1 to count($free)">
                            <xsl:variable name="outerIndex" select="."/>
                            <xsl:element name="step">
                                    <xsl:for-each select="1 to $atomCount">
                                        <xsl:variable name="innerIndex" select="."/>
                                        <xsl:variable name="forceIndex" select="(($outerIndex - 1) * $atomCount * 3) + (($innerIndex - 1) * 3)"/>                                       
                                        <xsl:element name="force">
                                            <xsl:if test="count($f) &gt;= $forceIndex + 3">
                                                <xsl:variable name="forceValues" select="($f[$forceIndex + 1], $f[$forceIndex + 2], $f[$forceIndex + 3])"/>
                                                <xsl:value-of select="helper:calcAxisLength($forceValues)" />
                                            </xsl:if>
                                        </xsl:element>
                                    </xsl:for-each>
                            </xsl:element>
                        </xsl:for-each>
                    </xsl:variable>
                                                                                                                                
 

                     <xsl-text>gibbsFreeEnergyList:</xsl-text>
                     <xsl:for-each select="$free">
                         <xsl:value-of select="."/>
                         <xsl-text> </xsl-text>
                     </xsl:for-each>
                     <xsl-text>#;#</xsl-text>
                    <xsl-text>deltaList:</xsl-text>
                    <xsl:for-each select="$delta">
                        <xsl:value-of select="."/>
                        <xsl-text> </xsl-text>
                    </xsl:for-each>
                    <xsl-text>#;#</xsl-text>
                    <xsl-text>forcesList:</xsl-text>
                    <xsl:for-each select="1 to count($free)">
                        <xsl:variable name="outerIndex" select="."/>
                        <xsl:choose>
                            <xsl:when test="empty($maxforces/*:step[$outerIndex]/*:force[1]/text())">' '</xsl:when>
                            <xsl:otherwise>
                                <xsl:value-of select="max($maxforces/*:step[$outerIndex]/*:force)" />
                                <xsl-text> </xsl-text>
                            </xsl:otherwise>
                        </xsl:choose>
                    </xsl:for-each>
                    <xsl-text>noentropyEnergyList:</xsl-text>
                    <xsl:for-each select="$noentropy">
                    <xsl:value-of select="."/>
                    <xsl-text> </xsl-text>
                    </xsl:for-each>
                    <xsl-text>#;#</xsl-text>
                    <xsl-text>e0List:</xsl-text>
                    <xsl:for-each select="$e0">
                        <xsl:value-of select="."/>
                        <xsl-text> </xsl-text>
                    </xsl:for-each>
                    <xsl-text>#;#</xsl-text>
                    <xsl-text>fermiList:</xsl-text>
                    <xsl:for-each select="1 to count($free)">
                        <xsl:variable name="outerIndex" select="."/>
                        <xsl:choose>
                            <xsl:when test="not(exists($fermi[$outerIndex]))">' '</xsl:when>
                            <xsl:otherwise>
                                <xsl:value-of select="$fermi[$outerIndex]"/>
                                <xsl-text> </xsl-text>
                            </xsl:otherwise>
                        </xsl:choose>
                    </xsl:for-each>                                                                  
                </xsl:otherwise>
            </xsl:choose>         
        </xsl:if>
        <xsl:if test="count($energies) > 1 ">
            <!-- Multiple jobs, we're dealing with NEB calculations -->
 
            <xsl:variable name="baseEnergy" select="number(tokenize($energies[1]/*[@dictRef='cc:e0'],'\s+')[last()])"/>                
            <xsl-text>nebEnergiesList:</xsl-text>
            <xsl:for-each select="$energies">
                <xsl:variable name="step" select="format-number(position()-1, '00')"/>
                <xsl:variable name="rel_e0" select="format-number(round-half-to-even(number(tokenize(./*[@dictRef='cc:e0'], '\s+')[last()]) - $baseEnergy,10),'0.0000')"/>
                <xsl:value-of select="concat($rel_e0,' ')"/>
            </xsl:for-each>
            <xsl-text>#;#</xsl-text>

            <xsl:value-of select="concat('nebEnergiesListUnits:',helper:printUnitSymbol($energies[1]/*[@dictRef='cc:e0']/@units),'#;#')"/>  
        </xsl:if>
        
    </xsl:template>
   
    <!-- Eigenvalues -->
 <!--   <xsl:template name="eigenvalues">
        <xsl:variable name="eigenvalues" select="//module[@dictRef='cc:calculation']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='eigenvalues']"/>        
        <xsl:if test="count($eigenvalues) = 1">
            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#eigenvalues-{generate-id($eigenvalues)}" style="cursor: pointer; cursor: hand;" >
                    <h4 class="panel-title">                        
                        Eigenvalues                           
                    </h4>
                </div>
                <div id="eigenvalues-{generate-id($eigenvalues)}" class="panel-collapse collapse">
                    <div class="panel-body">    
                        <div class="row bottom-buffer">
                            <div class="col-md-12">
                                <div class="row">                                
                                    <xsl:for-each select="$eigenvalues/list">
                                        <div class="col-lg-6 col-md-6 col-sm-12">
                                            <xsl:variable name="nodeId" select="generate-id(.)"/>
                                            <xsl:variable name="spin" select="./scalar[@dictRef='cc:spin']"></xsl:variable>
                                            <h4>Spin <xsl:value-of select="if(compare($spin,'1') = 0) then
                                                'alpha'
                                                else
                                                'beta'"/></h4>
                                            <div class="button-group form-inline mb-2">
                                                Kpoint <select class="form-control form-control-sm ml-2" id="kpointEigen-{generate-id(.)}" onchange="buildEigenOccupTable('eigenvaluesspin-{$nodeId}','{$spin}', this.value)">
                                                            <xsl:for-each select="tokenize(./array[@dictRef='cc:serial'],'\s+')">
                                                                <option value="{.}"><xsl:value-of select="."/></option>
                                                            </xsl:for-each>
                                                       </select>
                                            </div>
                                            <script type="text/javascript">
                                                <xsl:variable name="numberOfRows" select=" number(./array[@dictRef='cc:eigen']/@size) idiv number(./array[@dictRef='cc:serial']/@size)"/>
                                                <xsl:variable name="eigen" select="tokenize(./array[@dictRef='cc:eigen'],'\s+')"/>
                                                <xsl:variable name="occup" select="tokenize(./array[@dictRef='cc:occup'],'\s+')"/>
                                                <xsl:for-each select="tokenize(./array[@dictRef='cc:serial'],'\s+')">
                                                    <xsl:variable name="outerIndex" select="position()"/>                                       
                                                    var eigenValuesSpin<xsl:value-of select="$spin"/>Serial<xsl:value-of select="."/> = [<xsl:for-each select="1 to $numberOfRows"><xsl:variable name="innerIndex" select="."/>['<xsl:value-of select="$eigen[ ($outerIndex - 1) * $numberOfRows + $innerIndex ]"/>','<xsl:value-of select="$occup[ ($outerIndex - 1) * $numberOfRows + $innerIndex ]"/>']<xsl:if test="$innerIndex &lt; $numberOfRows">,</xsl:if></xsl:for-each>];
                                                </xsl:for-each>
                                                
                                                $(document).ready(function() {
                                                    $('#<xsl:value-of select="concat('kpointEigen-',generate-id(.))"/>').prop('selectedIndex', -1);                                                    
                                                });
                                            </script>
                                            <table class="display" id="eigenvaluesspin-{$nodeId}"></table>
                                        </div>                                    
                                    </xsl:for-each>
                                </div>
                                <script type="text/javascript">
                                    function buildEigenOccupTable(tableId,spin,kpoint){                                                       
                                        $('table#' + tableId).dataTable({
                                            "aaData":  window['eigenValuesSpin' + spin + 'Serial' + kpoint],
                                            "aoColumns": [ { "sTitle": "Eigenvalues"},{"sTitle": "Occupation"}], 
                                            "bFilter": false,
                                            "bPaginate": true,
                                            "bSort": false,
                                            "bInfo": true,
                                            "bDestroy": true
                                        });
                                                                              
                                        var jsTable = "javascript:showDownloadOptions('" + tableId + "');"
                                        $('div#downloadTable' + tableId).remove();
                                        $('<div id="downloadTable'+ tableId + '" class="text-right"><a class="text-right" href="' + jsTable +'"><span class="text-right fa fa-download"/></a></div>').insertBefore('div#' + tableId +'_wrapper'); 
                                     
                                    }
                                </script>
                            </div>
                        </div>
                    </div>
                </div>
            </div>    
        </xsl:if>              
    </xsl:template>-->
   

    <!-- Magnetization -->
    <xsl:template name="magnetization">
        <xsl:if test="exists($magnetizations)">   
            
            <xsl:for-each select="//module[@cmlx:templateRef='job']">  <!-- There can be uploaded calculations without -->
                <xsl:variable name="index" select="position()"/>
                <xsl:variable name="magnetization" select=".//module[@dictRef='cc:finalization']/propertyList/property/module[@cmlx:templateRef='magnetization']"/>                                        
                <xsl:choose>                                            
                    <xsl:when test="exists($magnetization)">  
                        <xsl:value-of select="concat('magnetization',$index,':')"/>
                        <xsl-text>serial coeff-s coeff-p coeff-d coeff-f total &#xa;</xsl-text>
                        <xsl:variable name="hasF" select="exists($magnetization/array[@dictRef='v:coefff'])"/>
                        <xsl:variable name="serial" select="tokenize($magnetization/array[@dictRef='cc:serial'],'\s+')"/>
                        <xsl:variable name="coeffs" select="tokenize($magnetization/array[@dictRef='v:coeffs'],'\s+')"/>
                        <xsl:variable name="coeffp" select="tokenize($magnetization/array[@dictRef='v:coeffp'],'\s+')"/>
                        <xsl:variable name="coeffd" select="tokenize($magnetization/array[@dictRef='v:coeffd'],'\s+')"/>
                        <xsl:variable name="coefff" select="tokenize($magnetization/array[@dictRef='v:coefff'],'\s+')"/>
                        <xsl:variable name="coefftotal" select="tokenize($magnetization/array[@dictRef='v:coefftotal'],'\s+')"/>           
                        <xsl:variable name="totals" select="$magnetization/list[@cmlx:templateRef='totals']"/>

                        <xsl:for-each select="1 to count($serial)">
                            <xsl:variable name="outerIndex" select="."/>
                            <xsl:value-of select="concat($serial[$outerIndex],' ')"/>
                            <xsl:value-of select="concat($coeffs[$outerIndex],' ')"/>
                            <xsl:value-of select="concat($coeffp[$outerIndex],' ')"/>
                            <xsl:value-of select="concat($coeffd[$outerIndex],' ')"/>
                            <xsl:choose>
                                <xsl:when test="$hasF">
                                    <xsl:value-of select="concat($coefff[$outerIndex],' ')"/>
                                </xsl:when>
                                <xsl:otherwise><xsl-text>None </xsl-text></xsl:otherwise>
                            </xsl:choose>
                    
                            <xsl:value-of select="concat($coefftotal[$outerIndex],'&#xa;')"/>
                        </xsl:for-each>
                    </xsl:when>
                </xsl:choose>
            </xsl:for-each>                                                      
        </xsl:if>       
    </xsl:template>    
 
    <!-- Vibrational frequencies -->
    <xsl:template name="vibrations">        
        <xsl:if test="exists($vibrations)">
            <xsl:variable name="frequency" select="tokenize($vibrations/array[@dictRef='cc:frequency'],'\s+')"/>
            <xsl:variable name="type" select="tokenize($vibrations/array[@dictRef='v:freqtype'],'\s+')"/>
            <xsl:variable name="cell">
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc1']"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc2']"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc3']"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc4']"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc5']"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$finalMolecules[1]/crystal/scalar[@id='sc6']"/>
            </xsl:variable>    
            <xsl:value-of select="concat('cellParsFreqs:',$cell,'#;#')"/>
            
            <xsl-text>frequencies: </xsl-text>
            <xsl:for-each select="1 to count($frequency)">
                <xsl:variable name="outerIndex" select="."/>
                <xsl:value-of select="concat($frequency[$outerIndex],' ')"/>
            </xsl:for-each>
            <xsl-text>#;#</xsl-text> 
        </xsl:if>       
    </xsl:template>
    

 

    <xsl:function name="vasp:getDivision">
       <xsl:param name="klabels" as="xs:string*" />
       <xsl:param name="index" as="xs:integer" />
       <xsl:param name="kdiv" as="xs:integer" />
       <xsl:choose>
           <xsl:when test="$index = 1"> <!-- First and last elements -->
               <xsl:element name="division">
                   <xsl:attribute name="value" select="0" />
                   <xsl:attribute name="label" select="$klabels[$index]" />
               </xsl:element>
           </xsl:when>
           <xsl:when test="$index = count($klabels)">
               <xsl:element name="division">
                   <xsl:attribute name="value" select="$kdiv" />
                   <xsl:attribute name="label" select="$klabels[$index]" />
               </xsl:element>
           </xsl:when>
           <xsl:otherwise>
               <xsl:if test="compare($klabels[$index], $klabels[$index + 1]) != 0">
                <xsl:choose>                           
                   <xsl:when test="compare($klabels[$index], $klabels[$index - 1]) = 0">
                       <xsl:element name="division">
                           <xsl:attribute name="value" select="$kdiv" />
                            <xsl:attribute name="label" select="$klabels[$index]"/>
                       </xsl:element>
                   </xsl:when>
                   <xsl:otherwise>
                       <xsl:if test="not(vasp:isDivisionStop($klabels, $index - 1))">
                           <xsl:element name="division">
                               <xsl:attribute name="value" select="$kdiv" />
                               <xsl:attribute name="label" select="concat($klabels[$index], '|', $klabels[$index +1])"/>
                               <xsl:attribute name="stop" select="'true'" />
                           </xsl:element>
                       </xsl:if>                    
                   </xsl:otherwise>
                </xsl:choose>
               </xsl:if>
           </xsl:otherwise>
       </xsl:choose>       
   </xsl:function>

    <xsl:function name="vasp:isDivisionStop" as="xs:boolean">
        <xsl:param name="klabels" />
        <xsl:param name="index" />       
        <xsl:value-of  select="compare($klabels[$index], $klabels[$index + 1]) != 0 and compare($klabels[$index], $klabels[$index - 1]) != 0 " />                           
    </xsl:function>
    

    
    <!-- Override default templates -->
    <xsl:template match="text()"/>
    <xsl:template match="*"/>
    
</xsl:stylesheet>
