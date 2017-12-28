	.file	"optkern.cpp"
	.intel_syntax noprefix
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"%lu %d \n"
	.text
	.p2align 4,,15
	.globl	_Z15accumtest4_7_10PdS_S_S_S_Pli
	.type	_Z15accumtest4_7_10PdS_S_S_S_Pli, @function
_Z15accumtest4_7_10PdS_S_S_S_Pli:
.LFB1082:
	.cfi_startproc
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	mov	rbp, rsp
	.cfi_def_cfa_register 6
	push	r15
	push	r14
	push	r13
	push	r12
	push	rbx
	and	rsp, -32
	add	rsp, 8
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	mov	eax, DWORD PTR [rbp+16]
	mov	QWORD PTR [rsp-24], rsi
	mov	QWORD PTR [rsp-32], rdx
	mov	QWORD PTR [rsp-40], rcx
	mov	QWORD PTR [rsp-48], r8
	lea	ebx, [rax-4096]
	mov	DWORD PTR [rsp-52], eax
	test	ebx, ebx
	mov	DWORD PTR [rsp-16], ebx
	jle	.L11
	mov	ecx, DWORD PTR [rsp-16]
	xor	r15d, r15d
	sub	eax, 4097
	xor	r11d, r11d
	mov	DWORD PTR [rsp-12], eax
	mov	r9, rsi
	vmovapd	ymm3, YMMWORD PTR .LC0[rip]
	vxorpd	xmm2, xmm2, xmm2
	sub	ecx, r15d
	test	ecx, ecx
	jle	.L12
.L19:
	mov	esi, DWORD PTR [rsp-12]
	mov	r14, QWORD PTR [rsp-24]
	mov	r13, QWORD PTR [rsp-32]
	mov	r12, QWORD PTR [rsp-40]
	mov	rbx, QWORD PTR [rsp-48]
	shr	esi, 2
	sub	r14, r11
	sal	rsi, 5
	sub	r13, r11
	sub	r12, r11
	lea	rsi, [r11+32+rsi]
	sub	rbx, r11
	.p2align 4,,10
	.p2align 3
.L9:
	lea	rdi, [r14+r11]
	lea	rcx, [r13+0+r11]
	lea	rdx, [r12+r11]
	lea	r10, [rbx+r11]
	xor	r8d, r8d
	.p2align 4,,10
	.p2align 3
.L10:
	vxorpd	xmm0, xmm0, xmm0
	xor	eax, eax
.L6:
	vmovapd	ymm1, YMMWORD PTR [rcx+rax]
	vmovapd	ymm9, YMMWORD PTR [rax+32+rcx]
	vsubpd	ymm5, ymm1, YMMWORD PTR [r9+rax]
	vbroadcastsd	ymm4, QWORD PTR [rax+rdi]
	vsubpd	ymm10, ymm9, YMMWORD PTR [rax+32+r9]
	vmovapd	ymm14, YMMWORD PTR [rax+64+rcx]
	vbroadcastsd	ymm8, QWORD PTR [rax+32+rdi]
	vsubpd	ymm15, ymm14, YMMWORD PTR [rax+64+r9]
	vbroadcastsd	ymm13, QWORD PTR [rax+64+rdi]
	vandpd	ymm6, ymm5, ymm3
	vmovapd	ymm5, YMMWORD PTR [rax+96+rcx]
	vandpd	ymm11, ymm10, ymm3
	vfmsub132pd	ymm6, ymm4, YMMWORD PTR [rdx+rax]
	vmaxpd	ymm7, ymm2, ymm6
	vbroadcastsd	ymm4, QWORD PTR [rax+96+rdi]
	sub	rax, -128
	vsubpd	ymm6, ymm5, YMMWORD PTR [rax-32+r9]
	vfmsub132pd	ymm11, ymm8, YMMWORD PTR [rax-96+rdx]
	vmaxpd	ymm12, ymm2, ymm11
	vfmadd132pd	ymm7, ymm0, ymm7
	vandpd	ymm0, ymm15, ymm3
	vfmsub132pd	ymm0, ymm13, YMMWORD PTR [rax-64+rdx]
	vmaxpd	ymm1, ymm2, ymm0
	vfmadd132pd	ymm12, ymm7, ymm12
	vandpd	ymm7, ymm6, ymm3
	vfmsub132pd	ymm7, ymm4, YMMWORD PTR [rax-32+rdx]
	vmaxpd	ymm0, ymm2, ymm7
	cmp	rax, 1024
	vfmadd132pd	ymm1, ymm12, ymm1
	vfmadd132pd	ymm0, ymm1, ymm0
	jne	.L6
	vmovapd	YMMWORD PTR [r10+r8], ymm0
	add	r8, 32
	cmp	r8, 256
	jne	.L10
	add	r11, 32
	cmp	r11, rsi
	jne	.L9
.L8:
	add	r15d, 32
	sub	DWORD PTR [rsp-12], 32
	add	r9, 256
	cmp	r15d, DWORD PTR [rsp-16]
	jge	.L16
	mov	ecx, DWORD PTR [rsp-16]
	mov	r11, rsi
	sub	ecx, r15d
	test	ecx, ecx
	jg	.L19
.L12:
	mov	rsi, r11
	jmp	.L8
.L16:
	vzeroupper
.L2:
	mov	edx, DWORD PTR [rsp-52]
	mov	edi, OFFSET FLAT:.LC1
	xor	eax, eax
	lea	rsp, [rbp-40]
	pop	rbx
	pop	r12
	pop	r13
	pop	r14
	pop	r15
	pop	rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	jmp	printf
.L11:
	.cfi_restore_state
	xor	esi, esi
	jmp	.L2
	.cfi_endproc
.LFE1082:
	.size	_Z15accumtest4_7_10PdS_S_S_S_Pli, .-_Z15accumtest4_7_10PdS_S_S_S_Pli
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC0:
	.long	4290772992
	.long	1105199103
	.long	4290772992
	.long	1105199103
	.long	4290772992
	.long	1105199103
	.long	4290772992
	.long	1105199103
	.ident	"GCC: (GNU) 4.8.5 20150623 (Red Hat 4.8.5-16)"
	.section	.note.GNU-stack,"",@progbits
