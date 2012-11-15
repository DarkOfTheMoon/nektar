<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fourier Single Mode Basis, P=7</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_SM.xml</parameters>
    <files>
        <file description="Session File">Test_SM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.99098e-14</value>
            <value variable="v" tolerance="1e-12">1.55055e-14</value>
            <value variable="w" tolerance="1e-12">2.56621e-14</value>
            <value variable="p" tolerance="1e-12">8.91393e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.23821e-13</value>
            <value variable="v" tolerance="1e-12">5.94171e-14</value>
            <value variable="w" tolerance="1e-12">1.28009e-13</value>
            <value variable="p" tolerance="1e-12">6.60807e-12</value>
        </metric>
    </metrics>
</test>