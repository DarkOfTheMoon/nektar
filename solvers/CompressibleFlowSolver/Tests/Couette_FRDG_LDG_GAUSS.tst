<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRDG advection and LDG diffusion, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0879609</value>
            <value variable="rhou" tolerance="1e-12">60.336</value>
            <value variable="rhov" tolerance="1e-8">0.227321</value>
            <value variable="E" tolerance="1e-12">4924.14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0736795</value>
            <value variable="rhou" tolerance="1e-12">61.0601</value>
            <value variable="rhov" tolerance="2e-6">0.262767</value>
            <value variable="E" tolerance="1e-12">4423.81</value>
        </metric>
    </metrics>
</test>


