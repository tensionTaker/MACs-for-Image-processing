module pg_bit(x, y, p, g);
	input x, y;
	output p, g;
	xor pblock (p, x, y);
	and gblock (g, x, y);
endmodule

module black_cell(g_ik, p_ik, g_k1j, p_k1j, Gij, Pij);
	input g_ik, p_ik, g_k1j, p_k1j;
	output Gij, Pij;
	wire w1;
	
	and g1 (w1, p_ik, g_k1j);
	or g2 (Gij, g_ik, w1);
	and g3 (Pij, p_ik, p_k1j);
endmodule

module gray_cell(g_ik, p_ik, g_k1j, Gij);
	input g_ik, p_ik, g_k1j;
	output Gij;
	
	wire w1;
	or g1 (Gij, g_ik, w1);
	and g2 (w1, p_ik, g_k1j);
endmodule


module KS_ADDER(clk, rst_n, a, b, c_in, sum, cout);
input clk, rst_n, c_in;
input [31:0]a;
input [31:0]b;

output reg [31:0]sum;
output reg cout;

reg [31:0]a_reg;  //  
reg [31:0]b_reg;  //    After 1st block
reg c_in_reg;     // 

wire [31:0]Yout; //Later after clk Yout become sum
wire cout_p; // this after clk becomes cout


// Bitwise PG generation

wire [32:0]bitwise_p;
wire [32:0]bitwise_g;

assign bitwise_p[0] = 0;
assign bitwise_g[0] = c_in_reg;

genvar i;
generate for(i = 1; i <= 32; i = i+1) begin : pg_bit_wise
	pg_bit pg_generation (a_reg[i-1], b_reg[i-1], bitwise_p[i], bitwise_g[i]);
	end
endgenerate


// Note bits are still numbered 0 to 31 but in generation of 
// PG bit wise and group PG we number the actual bit content 
// by 1 to 32 that that is p[1], g[1] denote 
// p and g bit wise of a[0], b[0]
// g[0], p[0] is for cin
// Bitwise PG generation ends here

// Layer1 of PG Propagation

wire [31:0]l1p;
wire [31:0]l1g;

genvar l1;
generate for(l1 = 2; l1 <= 31; l1 = l1 + 1) begin : l1b
		black_cell l1bcs (bitwise_g[l1], bitwise_p[l1], bitwise_g[l1-1], bitwise_p[l1-1], l1g[l1], l1p[l1]);
	end
endgenerate 

assign l1p[1] = 1'b0;
assign l1p[0] = 1'b0;

gray_cell l1gray1 (bitwise_g[1], bitwise_p[1], bitwise_g[0], l1g[1]);

buf l1buffer1 (l1g[0], bitwise_g[0]);

// Layer 1 ends here

// Layer 2 starts
//module black_cell(g_ik, p_ik, g_k1j, p_k1j, Gij, Pij);
//module gray_cell(g_ik, p_ik, g_k1j, Gij);
wire [31:0]l2p;
wire [31:0]l2g;

//l2 black cells
genvar l2;
generate for(l2 = 4; l2 <= 31; l2 = l2 + 1) begin : l2b
	black_cell l2bcs (l1g[l2], l1p[l2], l1g[l2-2], l1p[l2-2], l2g[l2], l2p[l2]);
	end
endgenerate

gray_cell l2gray1 (l1g[3], l1p[3], l1g[1], l2g[3]);
gray_cell l2gray2 (l1g[2], l1p[2], l1g[0], l2g[2]);
assign l2p[3] = 1'b0;
assign l2p[2] = 1'b0;
assign l2p[1] = 1'b0;

assign l2p[0] = 1'b0;
assign l2g[0] = l1g[0];

buf l2buffer (l2g[1], l1g[1]);
// Layer 2 ends here

// Layer 3
//module black_cell(g_ik, p_ik, g_k1j, p_k1j, Gij, Pij);
//module gray_cell(g_ik, p_ik, g_k1j, Gij);
wire [31:0]l3g;
wire [31:0]l3p;

genvar l3;
generate for(l3 = 8; l3 <= 31; l3 = l3 + 1) begin : l3b
		black_cell l3bcs (l2g[l3], l2p[l3], l2g[l3-4], l2p[l3-4], l3g[l3], l3p[l3]);
	end
endgenerate 

gray_cell l3gray1 (l2g[7], l2p[7], l2g[3], l3g[7]);
gray_cell l3gray2 (l2g[6], l2p[6], l2g[2], l3g[6]);
gray_cell l3gray3 (l2g[5], l2p[5], l2g[1], l3g[5]);
gray_cell l3gray4 (l2g[4], l2p[4], l2g[0], l3g[4]);

assign l3p[7] = 1'b0;
assign l3p[6] = 1'b0;
assign l3p[5] = 1'b0;
assign l3p[4] = 1'b0;

buf l3buffer1 (l3g[3], l2g[3]);
buf l3buffer2 (l3g[2], l2g[2]);
assign l3p[3] = 1'b0;
assign l3p[2] = 1'b0;

assign l3g[1] = l2g[1];
assign l3p[1] = 1'b0;
assign l3g[0] = l2g[0];
assign l3p[0] = 1'b0;

// Layer 3 ends here


//Layer 4 
wire [31:0]l4g;
wire [31:0]l4p;

genvar l4;
generate for(l4 = 16; l4 <= 31; l4 = l4 + 1) begin : l4b
		black_cell l4bcs (l3g[l4], l3p[l4], l3g[l4-8], l3p[l4-8], l4g[l4], l4p[l4]);
	end
endgenerate

genvar l4gc;
generate for(l4gc = 8; l4gc <= 15; l4gc = l4gc + 1) begin : l4gcell
		gray_cell l4gcs (l3g[l4gc], l3p[l4gc], l3g[l4gc - 8], l4g[l4gc]);
	end
endgenerate

genvar l4bf;
generate for(l4bf = 4; l4bf <= 7; l4bf = l4bf + 1) begin : buffersl4
		buf l4buffer (l4g[l4bf], l3g[l4bf]);
	end
endgenerate
assign l4g[3] = l3g[3];
assign l4g[2] = l3g[2];
assign l4g[1] = l3g[1];
assign l4g[0] = l3g[0];

assign l4p[15:0] = 16'b0;

// Layer 4 ends here

// Layer 5
//module gray_cell(g_ik, p_ik, g_k1j, Gij);
wire [31:0]l5g;

genvar l5;
generate for(l5 = 16; l5 <= 31; l5 = l5 + 1) begin : gcs
		gray_cell l5gcs (l4g[l5], l4p[l5], l4g[l5-16], l5g[l5]);
	end
endgenerate

genvar l5bf;
generate for(l5bf = 8; l5bf <= 15; l5bf = l5bf + 1) begin : buffersl5
		buf l5buffer (l5g[l5bf], l4g[l5bf]);
	end
endgenerate

assign l5g[7] = l4g[7];
assign l5g[6] = l4g[6];
assign l5g[5] = l4g[5];
assign l5g[4] = l4g[4];
assign l5g[3] = l4g[3];
assign l5g[2] = l4g[2];
assign l5g[1] = l4g[1];
assign l5g[0] = l4g[0];

// Layer 5 ends here

// Final Summing layer
wire last1;
and coutG1 (last1, bitwise_p[32], l5g[31]);
or coutG2 (cout_p, bitwise_g[32], last1);

genvar s_out;
generate for(s_out = 0; s_out <=  31; s_out = s_out + 1) begin : summing 
	xor sumready (Yout[s_out], bitwise_p[s_out+1], l5g[s_out]);
	end
endgenerate


always @(posedge clk or negedge rst_n) begin
	if(!rst_n) begin
		a_reg <= 32'b0;
		b_reg <= 32'b0;
		c_in_reg <= 1'b0;
		
		sum <= 32'b0;
		cout <= 1'b0;
	end
	else begin
		a_reg <= a;
		b_reg <= b;
		c_in_reg <= c_in;
		
		// Output block
		sum <= Yout;
		cout <= cout_p;
	end
end
endmodule


module multipler(
input clk,
input rst,
input [7:0] A,
input [7:0] B,
output [15:0] Product
    );
wire [3:0]  A_low,A_high,B_low,B_high;
assign A_low = A[3:0];
assign A_high = A[7:4];
assign B_low = B[3:0];
assign B_high = B[7:4];

wire [7:0] Product1,Product2,Product3,Product4;
assign Product1 = A_low * B_low;
assign Product2 = A_high * B_low;
assign Product3 = A_low * B_high;
assign Product4 = A_high * B_high;

reg[7:0] reg_Product1,reg_Product2,reg_Product3,reg_Product4;
wire[8:0] Sum1 = reg_Product2 + reg_Product3;
wire[12:0] Sum1_shift = {Sum1,4'b0};
wire[15:0] Product1_shift,Product4_shift;
assign Product1_shift = {8'b0,reg_Product1};
assign Product4_shift = {reg_Product4,8'b0};

reg[12:0] reg_Sum1_shift;
reg[15:0] reg_Product1_shift,reg_Product4_shift;


wire[15:0] PP1,PP2,PP3;
assign PP1 =reg_Product1_shift;
assign PP2 ={3'b0,reg_Sum1_shift};
assign PP3 =reg_Product4_shift;
wire [15:0] result = PP1 + PP2 + PP3;
reg[15:0] reg_Product;
always @(posedge clk) begin
        if (rst) begin
            // Clear pipeline
            reg_Product1        <= 8'b0;
            reg_Product2        <= 8'b0;
            reg_Product3        <= 8'b0;
            reg_Product4        <= 8'b0;

            reg_Sum1_shift  <= 13'b0;
            reg_Product1_shift <= 16'b0;
            reg_Product4_shift    <= 16'b0;

            reg_Product    <= 16'b0;

        end else begin

            // Stage 0 → Stage 1
            reg_Product1 <= Product1;
            reg_Product2 <= Product2;
            reg_Product3 <= Product3;
            reg_Product4 <= Product4;

            // Stage 1 → Stage 2
            reg_Sum1_shift  <= Sum1_shift;
            reg_Product1_shift <= Product1_shift;
            reg_Product4_shift    <= Product4_shift;

            // Stage 2 → Stage 3
            reg_Product <= result;

        end
    end
    assign Product =reg_Product;
endmodule


module MAC(
    input  [7:0]  in1,
    input  [7:0]  h,
    input  [31:0] in2,
    input         rst_n,
    input         clk,
    output [31:0] z_out
);

    wire [15:0] mul_out;

    multipler u_mul (
        .clk(clk),
        .rst(~rst_n),   // note: multiplier uses active-high rst
        .A(in1),
        .B(h),
        .Product(mul_out)
    );

   
    // Sign-extend multiplier output to 32 bits

    wire [31:0] mul_ext = {16'b0, mul_out};


    // Kogge-Stone ADDER
    // z_out = mul_ext + in2

    wire [31:0] sum_wire;
    wire        cout_wire;

    KS_ADDER u_adder (
        .clk(clk),
        .rst_n(rst_n),
        .a(mul_ext),
        .b(in2),
        .c_in(1'b0),
        .sum(sum_wire),
        .cout(cout_wire)
    );
    
    // Final registered output

    assign z_out = sum_wire;

endmodule

module dff_c #(
parameter integer N = 8
) (
input clk,
input rst_n, // active-low synchronous reset
input [N-1:0] d,
output reg [N-1:0] q
);


always @(posedge clk) begin
	if (!rst_n)	q <= {N{1'b0}};
	else	q <= d;
end
endmodule

module modified_MAC(
    input  [7:0]  in1,
    input  [7:0]  h,
    input  [31:0] in2,
    input         rst_n,
    input         clk,
    input         slow_clk,
    output [31:0] z_out
);
	wire [7:0]in1_wire;
	wire [31:0]z_out_wire;
	
	dff_c #(.N(8)) d1 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(in1),
		.q(in1_wire)
	);
	
	MAC original (
		.in1(in1_wire),
		.h(h),
		.in2(in2),
		.rst_n(rst_n),
		.clk(clk),
		.z_out(z_out_wire)
	);

	dff_c #(.N(32)) d2 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(z_out_wire),
		.q(z_out)
	);
endmodule

module clk_div6(
    input  clk,
    input  rst_n,
    output reg slow_clk
);
    reg [2:0] cnt;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt      <= 3'd0;
            slow_clk <= 1'b0;
        end
        else begin
            if (cnt == 3'd2) begin
                cnt      <= 3'd0;
                slow_clk <= ~slow_clk;   // Toggle every 3 cycles → final period = 6 cycles
            end
            else begin
                cnt <= cnt + 1'b1;
            end
        end
    end
endmodule
// Top Module
module Kernel(
    input  [7:0] pixel_in,
    input        clk,
    input        rst_n,
    input 	 rst_clk,
    output [31:0] z_out,
    output        slow_clk_out     // <-- for testbench
);

    // --------- Slow clock generator ----------
    wire slow_clk;

    clk_div6 u_div6 (
        .clk(clk),
        .rst_n(rst_clk),
        .slow_clk(slow_clk)
    );

    assign slow_clk_out = slow_clk;
reg [7:0]h_in[0:8];
// Kernel Weights
always@(rst_n) begin
	h_in[0] = 8'd2;
	h_in[1] = 8'd1;
	h_in[2] = 8'd2;
	
	h_in[3] = 8'd1;
	h_in[4] = 8'd2;
	h_in[5] = 8'd1;
	
	h_in[6] = 8'd2;
	h_in[7] = 8'd1;
	h_in[8] = 8'd2;
end

wire [31:0]mid_val[0:9];

modified_MAC M1 (
	.in1(pixel_in),
	.h(h_in[0]),
	.in2(32'b0),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[0])
);

genvar m_gen;
generate for(m_gen = 1; m_gen <= 2; m_gen = m_gen + 1) begin : Munits
	 	modified_MAC Mxs (
			.in1(pixel_in),
			.h(h_in[m_gen]),
			.in2(mid_val[m_gen - 1]),
			.rst_n(rst_n),
			.clk(clk),
			.slow_clk(slow_clk),
			.z_out(mid_val[m_gen])
		);
	end
endgenerate
	wire [31:0]l1b[0:2]; 
	
	// Layer 1 buffer
	dff_c #(.N(32)) b1l1 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(mid_val[2]),
		.q(l1b[0])
	);
	dff_c #(.N(32)) b2l1 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(l1b[0]),
		.q(l1b[1])
	);
	dff_c #(.N(32)) b3l1 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(l1b[1]),
		.q(l1b[2])
	);	

modified_MAC M4 (
	.in1(pixel_in),
	.h(h_in[3]),
	.in2(l1b[2]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[3])
);

modified_MAC M5 (
	.in1(pixel_in),
	.h(h_in[4]),
	.in2(mid_val[3]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[4])
);

modified_MAC M6 (
	.in1(pixel_in),
	.h(h_in[5]),
	.in2(mid_val[4]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[5])
);
	wire [31:0]l2b[0:2];

	dff_c #(.N(32)) b1l2 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(mid_val[5]),
		.q(l2b[0])
	);
	dff_c #(.N(32)) b2l2 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(l2b[0]),
		.q(l2b[1])
	);
	dff_c #(.N(32)) b3l2 (
		.clk(slow_clk),
		.rst_n(rst_n), 
		.d(l2b[1]),
		.q(l2b[2])
	);
modified_MAC M7 (
	.in1(pixel_in),
	.h(h_in[6]),
	.in2(l2b[2]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[6])
);

modified_MAC M8 (
	.in1(pixel_in),
	.h(h_in[7]),
	.in2(mid_val[6]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[7])
);

modified_MAC M9 (
	.in1(pixel_in),
	.h(h_in[8]),
	.in2(mid_val[7]),
	.rst_n(rst_n),
	.clk(clk),
	.slow_clk(slow_clk),
	.z_out(mid_val[8])
);
assign z_out = mid_val[8];
endmodule
