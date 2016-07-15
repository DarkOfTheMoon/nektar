<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P, periodic BC</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=DirectMultiLevelStaticCond KovaFlow_varP_per.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_varP_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.10712e-05</value>
            <value variable="v" tolerance="1e-12">0.000107148</value>
	    <value variable="p" tolerance="1e-12">0.000558057</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00010879</value>
            <value variable="v" tolerance="1e-12">0.000368186</value>
	    <value variable="p" tolerance="1e-12">0.00334077</value>
        </metric>
    </metrics>
</test>

