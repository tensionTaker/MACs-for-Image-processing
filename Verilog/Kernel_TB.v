`timescale 1ns/1ps

module Kernel_TB;

    reg clk = 0;


    // Clock periods
    parameter FAST_PERIOD = 10;     // fast clock = 100 MHz
    parameter SLOW_PERIOD = 60;     // slow clock = fast/6

    // generate clocks
    always #(FAST_PERIOD/2) clk = ~clk;
    //always #(SLOW_PERIOD        repeat (80) @(posedge clk);/2) slow_clk = ~slow_clk;

    reg rst_n;
    reg rst_clk;
    reg [7:0] pixel_in;
    wire [31:0] z_out;
    wire slow_clk;
    // DUT instance
    Kernel DUT (
        .pixel_in(pixel_in),
        .clk(clk),
        .rst_n(rst_n),
        .rst_clk(rst_clk),
        .z_out(z_out), 
        .slow_clk_out(slow_clk)
    );

    integer i;

    initial begin
        // GTKWave VCD dump
        $dumpfile("tb_kernel.vcd");
        $dumpvars(0, Kernel_TB);

        // initial values
        rst_n = 0;
        rst_clk = 0;
        pixel_in = 0;
	repeat (3) @(posedge clk);
	rst_clk = 1;
        // hold reset for a few fast cycles
        repeat (80) @(posedge clk);
        rst_n = 1;

        // apply a long sequence of pixels on slow clock
        for (i = 0; i < 200; i = i + 1) begin
            @(posedge slow_clk);
            #30;
            pixel_in = i;
        end

        // let simulation continue a bit for pipeline to flush
        repeat (200) @(posedge clk);

        $finish;
    end

endmodule
