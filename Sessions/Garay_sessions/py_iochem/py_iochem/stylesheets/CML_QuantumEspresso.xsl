

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
		xmlns:qex="http://www.iochem-bd.org/dictionary/qespresso/"
		xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"

		xpath-default-namespace="http://www.xml-cml.org/schema"
		exclude-result-prefixes="xs xd cml ckbk qex helper cmlx"
		version="2.0">

		<xd:doc scope="stylesheet">
				<xd:desc>
						<xd:p><xd:b>Created on:</xd:b>May 29, 2024</xd:p>
						<xd:p><xd:b>Author:</xd:b>Moisés Álvarez Moreno, Diego Garay Ruiz</xd:p>
						<xd:p><xd:b>Center:</xd:b>Institute of Chemical Research of Catalonia</xd:p>
						<xd:p><xd:b>Center:</xd:b>Universitat Rovira i Virgili</xd:p>
				</xd:desc>
		</xd:doc>
		<xsl:include href="chem_helper.xsl"/>
		<xsl:include href="qespresso_helper.xsl"/>
		<xsl:output method="text" omit-xml-declaration="yes" indent="no" />
		<xsl:strip-space elements="*"/>

		<xsl:param name="title"/>
		<xsl:param name="author"/>

		<!-- Environment module -->


			<xsl:variable name="programParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:program']"/>
			<xsl:variable name="versionParameter" select="//module[@id='job'][1]/module[@id='environment']/parameterList/parameter[@dictRef='cc:programVersion']"/>

		<!-- Initialization module -->
			<xsl:variable name="modNames" select="//cml:module[@dictRef='cc:job']/module[@id='environment']/parameterList/parameter[@dictRef='cc:module']/scalar/text()"/>
			<xsl:variable name="calcParameter" select="replace(qex:getInputParameter(., 'CONTROL', 'calculation'), '&quot;','')"/>
			<xsl:variable name="calcType" select="qex:getCalcType($modNames,$calcParameter)"/>

			<xsl:variable name="otherComponents" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']/module[@id='otherComponents']"/>
			<xsl:variable name="userDefinedModules" select="$otherComponents/module[@dictRef='cc:userDefinedModule']"/>

			<xsl:variable name="species" select="tokenize(//cml:module[@cmlx:templateRef='species']/list[@cmlx:templateRef='species']/array[@dictRef='qex:specie'],'\s+')"/>
			<xsl:variable name="atomSpecies" select="tokenize((//module[@id='job'][1]/module[@dictRef='cc:initialization']/parameterList/parameter[@dictRef='qex:specie']/array)[1], '\s+')"/>
			<xsl:variable name="specieToPseudopotential">
					<xsl:for-each select="$species">
							<xsl:variable name="specie" select="."/>
							<xsl:variable name="outerIndex" select="position()"/>
							<xsl:element name="specie" namespace="http://www.xml-cml.org/schema">
									<xsl:attribute name="specie" select="$specie"/>
									<xsl:attribute name="basis" select="qex:getBasis($userDefinedModules/module[@cmlx:templateRef='pseudopotential' and child::scalar[@dictRef='cc:serial' and text() = $outerIndex]]/scalar[@dictRef='qex:pseudopotential'])"/>
									<xsl:attribute name="pseudofile" select="$userDefinedModules/module[@cmlx:templateRef='pseudopotential' and child::scalar[@dictRef='cc:serial' and text() = $outerIndex]]/scalar[@dictRef='qex:pseudofile']"/>
							</xsl:element>
					</xsl:for-each>
			</xsl:variable>
			<xsl:variable name="pseudoFilenames" select="$userDefinedModules/module[@cmlx:templateRef='pseudopotential']/scalar[@dictRef='qex:pseudofile']/text()"/>

			<xsl:variable name="functionals" select="helper:trim(replace($otherComponents/module[@cmlx:templateRef='parameters']//scalar[@dictRef='cc:parameter' and text() = 'Exchange-correlation']/following-sibling::scalar[@dictRef='cc:value']/text(),'\(.*\)',''))"/>
			<xsl:variable name="functionalsFromFilenames" select="distinct-values(qex:getFunctionalsFromFilenames($pseudoFilenames))"/>
			<xsl:variable name="functionalMethods" select="distinct-values(qex:getMethodsFromFilenames($pseudoFilenames))"/>
			<xsl:variable name="basis" select="distinct-values(qex:getBasis($userDefinedModules/module[@cmlx:templateRef='pseudopotential']/scalar[@dictRef='qex:pseudopotential']/text()))"/>

			<xsl:variable name="kpoints" select="$otherComponents/module[@cmlx:templateRef='kpoints']"/>
			<xsl:variable name="ldau" select="$otherComponents/module[@cmlx:templateRef='ldau']"/>
			<xsl:variable name="parameters" select="$otherComponents/module[@cmlx:templateRef='parameters']/list"/>
			<xsl:variable name="pointGroup" select="$otherComponents/module[@cmlx:templateRef='point.group']/scalar[@dictRef='cc:pointgroup']" />

			<xsl:variable name="temperature" select="$otherComponents/module[@cmlx:templateRef='environ']//scalar[@dictRef='cc:parameter' and text() = 'static permittivity']/following-sibling::scalar[@dictRef='cc:value']"/>
			<xsl:variable name="pressure" select="$otherComponents/module[@cmlx:templateRef='environ']//scalar[@dictRef='cc:parameter' and text() = 'external pressure in input']/following-sibling::scalar[@dictRef='cc:value']"/>
		<!-- Finalization module -->
		<xsl:variable name="frequencies" select="//module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='frequencies']"/>
		<xsl:variable name="absorption" select="//module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='absorptionspec']"/>
		<xsl:variable name="bands" select="//module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='bands']" />
		<xsl:variable name="pdos" select="//module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='pdos']" />
		<xsl:variable name="phonon" select="//module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='phonon']" />
		<xsl:variable name="method" select="if(exists($absorption)) then 'TDDFT' else 'DFT'"/>


		<!-- Geometry -->
		<xsl:variable name="energies" select="//module[@id='job']/module[@dictRef='cc:finalization']/module[@dictRef='cc:userDefinedModule']/module[@cmlx:templateRef='energies']"/>
		<xsl:variable name="finalMolecule" select="//module[@dictRef='cc:finalization' and child::molecule]//molecule"/>

		<xsl:template match="/">
			<xsl:call-template name="generalInfo"/>
			<xsl:call-template name="settings"/>
			<xsl:call-template name="atomicCoordinates">
				<xsl:with-param name="molecule" select="$finalMolecule"/>
			</xsl:call-template>

			<xsl:call-template name="ldau"/>
			<xsl:call-template name="pointGroup"/>
			<xsl:call-template name="kpoints"/>

			<xsl:for-each select="//cml:module[@dictRef='cc:job']">
					<xsl:variable name="job" select="."/>
					<xsl:variable name="moduleName" select="$job/module[@id='environment']/parameterList/parameter[@dictRef='cc:module']"/>
					<xsl:variable name="finalOtherComponents" select="$job/module[@dictRef='cc:finalization']/module[@id='otherComponents']" />
					<xsl:call-template name="force"> <xsl:with-param name="forces" select="$finalOtherComponents/module[@cmlx:templateRef='forces']"/> </xsl:call-template>
					<xsl:call-template name="energy"> <xsl:with-param name="energy" select="$finalOtherComponents/module[@cmlx:templateRef='energies']"/> </xsl:call-template>
					<xsl:call-template name="magnetization"><xsl:with-param name="magnetizations" select="$finalOtherComponents/module[@cmlx:templateRef='magnetic']"/></xsl:call-template>
					<xsl:call-template name="frequencies"/> 
					<xsl:text> FETCH_JOB_END&#xa;</xsl:text>
				
			</xsl:for-each>
		</xsl:template>

		<!-- General Info -->
		<xsl:template name="generalInfo">

				<xsl:value-of select="concat('program:',$programParameter,'#;#')"/>
				<xsl:value-of select="concat('version:',$versionParameter,'#;#')"/>
				<xsl:if test="$author">
						<xsl:value-of select="concat('author:',$author,'#;#')"/>
				</xsl:if>
				<xsl:variable name="calcFormula" select="helper:calculateHillNotationFormula($finalMolecule/atomArray/atom)"/>
				<xsl:value-of select="concat('formula:',string-join($calcFormula),'#;#')"/>
				<xsl:value-of select="concat('calctype:',$calcType,'#;#')"/>
				<xsl:value-of select="concat('method:',$method,'#;#')"/>
				<xsl:variable name="selFunctional" select="if(exists($functionals)) then $functionals else $functionalsFromFilenames"/>
				<xsl:value-of select="concat('functional:',$selFunctional,'#;#')"/>

				<xsl:if test="exists($temperature)">
					<xsl:value-of select="concat('temperature:',$temperature,'#;#')"/>
				</xsl:if>
				<xsl:if test="exists($pressure)">
					<xsl:variable name="pressureRounded" select="round-half-to-even($pressure/text(),2)"/>
					<xsl:variable name="pressUnits" select="helper:printUnitSymbol($pressure/@units)"/>
					<xsl:value-of select="concat('pressure:',$pressureRounded,'#;#')"/>
					<xsl:value-of select="concat('pressUnits:',$pressUnits,'#;#')"/>
				</xsl:if>
		</xsl:template>

		<!-- Settings -->
		<xsl:template name="settings">
			<xsl:for-each select="$parameters">
				<xsl:variable name="parameter" select="."/>
					<xsl:value-of select="concat($parameter/scalar[@dictRef='cc:parameter'],':',$parameter/scalar[@dictRef='cc:value'],'#;#')"/>
					<xsl:value-of select="concat($parameter/scalar[@dictRef='cc:parameter'],'Units:',$parameter/scalar[@dictRef='cc:value']/@units,'#;#')"/>
			</xsl:for-each>
		</xsl:template>

		<!-- Atomic coordinates -->
		<xsl:template name="atomicCoordinates">
				<xsl:param name="molecule"/>
				<xsl:value-of select="concat('cell_a:',$molecule/crystal/scalar[@title='a'],'#;#')"/>
				<xsl:value-of select="concat('cell_b:',$molecule/crystal/scalar[@title='b'],'#;#')"/>
				<xsl:value-of select="concat('cell_c:',$molecule/crystal/scalar[@title='c'],'#;#')"/>
				<xsl:value-of select="concat('cell_alpha:',$molecule/crystal/scalar[@title='alpha'],'#;#')"/>
				<xsl:value-of select="concat('cell_beta:',$molecule/crystal/scalar[@title='beta'],'#;#')"/>
				<xsl:value-of select="concat('cell_gamma:',$molecule/crystal/scalar[@title='gamma'],'#;#')"/>
				<!-- Repeat the loop to have Cartesian and fractional coordinates as well as additional parameters -->
				<xsl-text> geometryCartesian: </xsl-text>
				<xsl:for-each select="$molecule/cml:atomArray/cml:atom">
						<xsl:variable name="elementType" select="./@elementType"/>
						<xsl:value-of select="concat($elementType,' ',@x3,' ',@y3,' ',@z3,'&#xa;')"/>
				</xsl:for-each>
				<xsl-text> #;# </xsl-text>
				<xsl-text> geometryFractional: </xsl-text>
				<xsl:for-each select="$molecule/cml:atomArray/cml:atom">
						<xsl:variable name="elementType" select="./@elementType"/>
						<xsl:value-of select="concat($elementType,' ',@xFract,' ',@yFract,' ',@zFract,'&#xa;')"/>
				</xsl:for-each>
				<xsl-text> #;# </xsl-text>
				<xsl-text> speciesPseudopotentials: </xsl-text>
				<xsl:for-each select="$molecule/cml:atomArray/cml:atom">
						<xsl:variable name="elementType" select="./@elementType"/>
						<xsl:variable name="id" select="@id"/>
						<xsl:variable name="outerIndex" select="position()"/>
						<xsl:variable name="specie" select="$atomSpecies[$outerIndex]"/>
						<xsl:variable name="basis" select="$specieToPseudopotential/specie[@specie = $atomSpecies[$outerIndex]]/@basis"/>
						<xsl:variable name="pseudo" select="$specieToPseudopotential/specie[@specie = $atomSpecies[$outerIndex]]/@pseudofile"/>
						<xsl:value-of select="concat($specie,' ',$basis,' ',$pseudo,'&#xa;')"/>
				</xsl:for-each>
				<xsl-text> #;# </xsl-text>
		</xsl:template>

		<!-- kpoint list -->
		<xsl:template name="kpoints">
				<xsl:if test="exists($kpoints)">
						<xsl:variable name="alat" select="//module[@id='job'][1]/module[@dictRef='cc:initialization']//module[@id='otherComponents']/module[@cmlx:templateRef='parameters']/list/scalar[@dictRef='cc:parameter' and text() = 'lattice parameter (alat)']/following-sibling::scalar[@dictRef='cc:value']"/>
						<xsl:variable name="twopi_alat2angstrom" select="(0.529177249) div (2 * 3.14159265359 * $alat)"/>
						<xsl:variable name="subdivisionN" select="tokenize($kpoints/array[@dictRef='qex:subdivisionN'],'\s+')"/>
						<xsl:variable name="shiftS" select="tokenize($kpoints/array[@dictRef='qex:shiftS'],'\s+')"/>

						<xsl:variable name="kpointCrystalNumber" select="$kpoints/matrix[@dictRef='qex:kpointlist']/@rows"/>
						<xsl:variable name="kpointCartesian" select="tokenize($kpoints/matrix[@dictRef='qex:kpointlist.cartesian'],'\s+')"/>

						<xsl:variable name="kpointCartesianNumber" select="$kpoints/matrix[@dictRef='qex:kpointlist.cartesian']/@rows"/>
						<xsl:variable name="kpointCrystal" select="tokenize($kpoints/matrix[@dictRef='qex:kpointlist'],'\s+')"/>

						<xsl:variable name="kpointWeight" select="tokenize($kpoints/array[@dictRef='qex:kpointweight'][last()],'\s+')"/>
						<xsl:value-of select="concat('kpointScheme:',$kpoints/scalar[@dictRef='qex:meshScheme'],'#;#')"/>
						<xsl-text>kpointSubdivisions: </xsl-text>
						<xsl:for-each select="1 to count($subdivisionN)"><xsl:variable name="outerIndex" select="."/><xsl:value-of select="$subdivisionN[$outerIndex]"/></xsl:for-each>
						<xsl-text>#;#</xsl-text>
						<xsl-text>kpointShifts: </xsl-text>
						<xsl:for-each select="1 to count($shiftS)"><xsl:variable name="outerIndex" select="."/><xsl:value-of select="$shiftS[$outerIndex]"/></xsl:for-each>
						<xsl-text>#;#</xsl-text>
						<xsl-text>kpointListCartesian:</xsl-text>
						<xsl:for-each select="1 to $kpointCartesianNumber">
								<xsl:variable name="outerIndex" select="."/>
										<xsl:variable name="kpointCartx" select="format-number($twopi_alat2angstrom * number($kpointCartesian[3 * ($outerIndex - 1) + 1]),'#0.0000000')"/>
										<xsl:variable name="kpointCarty" select="format-number($twopi_alat2angstrom * number($kpointCartesian[3 * ($outerIndex - 1) + 2]),'#0.0000000')"/>
										<xsl:variable name="kpointCartz" select="format-number($twopi_alat2angstrom * number($kpointCartesian[3 * ($outerIndex - 1) + 3]),'#0.0000000')"/>
										<xsl:value-of select="concat($kpointCartx,' ',$kpointCarty,' ',$kpointCartz,' ',$kpointWeight[$outerIndex],'&#xa;')"/>
						</xsl:for-each>
						<xsl-text>#;#</xsl-text>
						<xsl:if test="exists($kpointCrystal)">
							<xsl-text>kpointListCrystal:</xsl-text>
							<xsl:for-each select="1 to $kpointCrystalNumber">
							<xsl:variable name="outerIndex" select="."/>
									<xsl:variable name="kpointCrysx" select="$kpointCrystal[3 * ($outerIndex - 1) + 1]"/>
									<xsl:variable name="kpointCrysy" select="$kpointCrystal[3 * ($outerIndex - 1) + 2]"/>
									<xsl:variable name="kpointCrysz" select="$kpointCrystal[3 * ($outerIndex - 1) + 3]"/>
									<xsl:value-of select="concat($kpointCrysx,' ',$kpointCrysy,' ',$kpointCrysz,' ',$kpointWeight[$outerIndex],'&#xa;')"/>
							</xsl:for-each>
							<xsl-text>#;#</xsl-text>
						</xsl:if>
				</xsl:if>
		</xsl:template>

		<!-- LDA+U table -->
		<xsl:template name="ldau">
				<xsl:if test="exists($ldau)">
						<xsl:variable name="ldauSpecie" select="tokenize($ldau//array[@dictRef='qex:specie'],'\s+')"/>
						<xsl:variable name="ldauL" select="tokenize($ldau//array[@dictRef='qex:l'],'\s+')"/>
						<xsl:variable name="ldauU" select="tokenize($ldau//array[@dictRef='qex:u'],'\s+')"/>
						<xsl:variable name="ldauAlpha" select="tokenize($ldau//array[@dictRef='qex:alpha'],'\s+')"/>
						<xsl:variable name="ldauJ0" select="tokenize($ldau//array[@dictRef='qex:j0'],'\s+')"/>
						<xsl:variable name="ldauBeta" select="tokenize($ldau//array[@dictRef='qex:beta'],'\s+')"/>
						<xsl-text> ldauBlock: </xsl-text>
						<xsl:for-each select="1 to count($ldauSpecie)">
								<xsl:variable name="outerIndex" select="."/>
										<xsl:value-of select="concat($ldauSpecie[$outerIndex],' ')"/>
										<xsl:value-of select="concat($ldauL[$outerIndex],' ')"/>
										<xsl:value-of select="concat($ldauU[$outerIndex],' ')"/>
										<xsl:value-of select="concat($ldauAlpha[$outerIndex],' ')"/>
										<xsl:value-of select="concat($ldauJ0[$outerIndex],' ')"/>
										<xsl:value-of select="concat($ldauBeta[$outerIndex],'&#xa;')"/>
						</xsl:for-each>
						<xsl-text> #;# </xsl-text>
				</xsl:if>
		</xsl:template>

		<!-- Point group line -->
		<xsl:template name="pointGroup">
				<xsl:if test="exists($pointGroup)">
					<xsl:value-of select="concat('pointGroup:',$pointGroup,'#;#')"/>
				</xsl:if>
		</xsl:template>

		<!-- Forces section -->
		<xsl:template name="force">
				<xsl:param name="forces"/>
				<xsl:if test="exists($forces)">
					<xsl:variable name="force" select="$forces/scalar[@dictRef='cc:force']"/>
					<xsl:variable name="scfCorrection" select="$forces/scalar[@dictRef='qex:scfCorrection']"/>
					<xsl:variable name="dispEnergy" select="$forces/scalar[@dictRef='qex:dispenergy']"/>
					<xsl:if test="exists($force)">
						<xsl:value-of select="concat('totalForce:',$force,'#;#')"/>
						<xsl:value-of select="concat('totalForceUnits:',helper:printUnitSymbol($force/@units),'#;#')"/>
					</xsl:if>
					<xsl:if test="exists($scfCorrection)">
							<xsl:value-of select="concat('totalSCFCorrection:',$scfCorrection,'#;#')"/>
							<xsl:value-of select="concat('totalSCFCorrectionUnits:',helper:printUnitSymbol($scfCorrection/@units),'#;#')"/>
					</xsl:if>
					<xsl:if test="exists($dispEnergy)">
							<xsl:value-of select="concat('totalDispersionForce:',$dispEnergy,'#;#')"/>
							<xsl:value-of select="concat('totalDispersionForceUnits:',helper:printUnitSymbol($dispEnergy/@units),'#;#')"/>
					</xsl:if>
				</xsl:if>
		</xsl:template>

		<!-- Final energy -->
		<xsl:template name="energy">
				<xsl:param name="energy"/>
				<xsl:if test="exists($energy)">
						<xsl:variable name="energyLabels">
								<variable name="qex:fermiener" label="Fermi energy"/>
								<variable name="qex:totalener" label="Total energy"/>
								<variable name="qex:harrisfoulkes" label="Harris-Foulkes estimate"/>
								<variable name="qex:sscfaccuracy" label="Estimated scf accuracy"/>
								<variable name="qex:oneelec" label="One-electron contribution"/>
								<variable name="qex:hartee" label="Hartree contribution"/>
								<variable name="qex:xc" label="XC contribution"/>
								<variable name="qex:ewald" label="Ewald contribution"/>
								<variable name="qex:dispcorr" label="Dispersion Correction"/>
								<variable name="qex:hubbardener" label="Hubbard energy"/>
								<variable name="qex:smearing" label="Smearing contrib. (-TS)"/>
								<variable name="qex:totalmag" label="Total magnetization"/>
								<variable name="cc:absolutemag" label="Absolute magnetization"/>
						</xsl:variable>

						<xsl:for-each select="$energy/scalar">
								<xsl:variable name="energyField" select="."/>
								<xsl:variable name="energyFieldLabel" select="$energyLabels/*:variable[@name=$energyField/@dictRef]/@label"/>
								<xsl:variable name="energyFieldUnits" select="helper:printUnitSymbol($energyField/@units)"/>
								<xsl:value-of select="concat($energyFieldLabel,':',$energyField,'#;#')" />
								<xsl:value-of select="concat($energyFieldLabel,'Units:',$energyFieldUnits,'#;#')" />
						</xsl:for-each>
						<xsl:if test="$calcType = $qex:GeometryOptimization and (count(//module[@cmlx:templateRef='energies']) &gt; 1)">
							<xsl-text> finalEnergies: </xsl-text>
							<xsl:variable name="finalEnergies" select="//module[@cmlx:templateRef='energies']" />
								<xsl:for-each select="1 to count($finalEnergies)">
									<xsl:variable name="outerIndex" select="."/>
									<xsl:value-of select="concat($finalEnergies[$outerIndex]/scalar[@dictRef='qex:totalener'] * 0.073498617649,' ')"/>
									<xsl:if test="$outerIndex &lt; count($finalEnergies)">
									</xsl:if>
								</xsl:for-each>
								<xsl-text> #;# </xsl-text>
							</xsl:if>
				</xsl:if>
		</xsl:template>

		<!-- Magnetization -->
		<xsl:template name="magnetization">
				<xsl:param name="magnetizations"/>
				<xsl:if test="exists($magnetizations)">
						<xsl:variable name="serial" select="tokenize($magnetizations/array[@dictRef='cc:serial'],'\s+')"/>
						<xsl:variable name="charge" select="tokenize($magnetizations/array[@dictRef='qex:charge'],'\s+')"/>
						<xsl:variable name="magn" select="tokenize($magnetizations/array[@dictRef='qex:magn'],'\s+')"/>
						<xsl:variable name="constr" select="tokenize($magnetizations/array[@dictRef='qex:constr'],'\s+')"/>
						<xsl-text> magnetization: </xsl-text>
						<xsl-text> index atom charge magnet constr '&#xa;'</xsl-text>
						<xsl:for-each select="1 to count($serial)">
							<xsl:variable name="outerIndex" select="."/>
							<xsl:value-of select="concat($serial[$outerIndex],' ')"/>
							<xsl:value-of select="concat($finalMolecule/atomArray/atom[$outerIndex]/@elementType,' ')"/>
							<xsl:value-of select="concat($charge[$outerIndex],' ')"/>
							<xsl:value-of select="concat($magn[$outerIndex],' ')"/>
							<xsl:value-of select="concat($constr[$outerIndex],'&#xa;')"/>
							<!-- <xsl:if test="$outerIndex != count($serial)">,</xsl:if> -->
						</xsl:for-each>
					<xsl-text> #;# </xsl-text>
				</xsl:if>
		</xsl:template>



		<xsl:template name="frequencies">
				<xsl:param name="section"/>
				<xsl:if test="exists($frequencies)">
						<xsl:variable name="serial" select="tokenize($frequencies/array[@dictRef='cc:serial'],'\s+')"/>
						<xsl:variable name="values" select="tokenize($frequencies/array[@dictRef='cc:frequency'],'\s+')"/>
						<xsl-text>frequencies:</xsl-text>
						<xsl:for-each select="1 to count($serial)">
								<xsl:variable name="outerIndex" select="."/>
										<xsl:value-of select="concat($values[$outerIndex],' ')"/>
						</xsl:for-each>
					<xsl-text>#;#</xsl-text>
				</xsl:if>
		</xsl:template>



		<xsl:function name="qex:formatXaxisLabel">
				<xsl:param name="label" />

				<xsl:variable name="l" select="replace($label,'!', '')" />
				<xsl:choose>
						<xsl:when test="matches($l,'G')">Γ</xsl:when>
						<xsl:otherwise><xsl:value-of select="$l"/></xsl:otherwise>
				</xsl:choose>

		</xsl:function>



		<!-- Override default templates -->
		<xsl:template match="text()"/>
		<xsl:template match="*"/>

</xsl:stylesheet>
