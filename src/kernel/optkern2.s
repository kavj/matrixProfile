	.file	"optkern2.cpp"
	.intel_syntax noprefix
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"%lu %d \n"
	.text
	.p2align 4,,15
	.globl	_Z14accumtest4_7_9PdS_S_S_S_Pli
	.type	_Z14accumtest4_7_9PdS_S_S_S_Pli, @function
_Z14accumtest4_7_9PdS_S_S_S_Pli:
.LFB1081:
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
	sub	rsp, 16
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	mov	eax, DWORD PTR [rbp+16]
	mov	QWORD PTR [rsp-96], rsi
	mov	QWORD PTR [rsp-104], rdx
	mov	QWORD PTR [rsp-112], rcx
	mov	QWORD PTR [rsp-120], r8
	lea	ebx, [rax-4096]
	mov	DWORD PTR [rsp-88], eax
	test	ebx, ebx
	mov	DWORD PTR [rsp-84], ebx
	jle	.L8
	mov	edx, DWORD PTR [rsp-84]
	xor	r14d, r14d
	xor	r10d, r10d
	lea	r15d, [rax-4097]
	mov	r11, rsi
	vmovapd	ymm10, YMMWORD PTR .LC0[rip]
	vxorpd	xmm9, xmm9, xmm9
	sub	edx, r14d
	test	edx, edx
	jle	.L9
.L14:
	mov	ecx, r15d
	mov	r13, QWORD PTR [rsp-104]
	mov	r12, QWORD PTR [rsp-112]
	shr	ecx, 2
	mov	rbx, QWORD PTR [rsp-96]
	mov	r9, QWORD PTR [rsp-120]
	sal	rcx, 5
	lea	rsi, [r10+32+rcx]
	sub	r13, r10
	sub	r12, r10
	sub	rbx, r10
	.p2align 4,,10
	.p2align 3
.L7:
	vxorpd	xmm15, xmm15, xmm15
	lea	rcx, [r12+r10]
	lea	rax, [r13+0+r10]
	lea	r8, [rbx+r10]
	mov	rdx, r11
	xor	edi, edi
	vmovapd	ymm14, ymm15
	vmovapd	ymm13, ymm15
	vmovapd	YMMWORD PTR [rsp-80], ymm15
	vmovapd	ymm12, ymm15
	vmovapd	ymm11, ymm15
	vmovapd	YMMWORD PTR [rsp-48], ymm15
	vmovapd	YMMWORD PTR [rsp-16], ymm15
	.p2align 4,,10
	.p2align 3
.L6:
	vmovapd	ymm8, YMMWORD PTR [rax]
	add	rdx, 32
	add	rax, 32
	vmovapd	ymm5, YMMWORD PTR [rax+32]
	add	rcx, 32
	vsubpd	ymm2, ymm8, YMMWORD PTR [rdx-32]
	vmovapd	ymm7, YMMWORD PTR [rax]
	vsubpd	ymm1, ymm5, YMMWORD PTR [rdx+32]
	vbroadcastsd	ymm0, QWORD PTR [r8+rdi]
	add	rdi, 32
	vsubpd	ymm6, ymm7, YMMWORD PTR [rdx]
	vandpd	ymm3, ymm2, ymm10
	vandpd	ymm2, ymm1, ymm10
	vmovapd	ymm1, YMMWORD PTR [rax+96]
	vfmsub132pd	ymm3, ymm0, YMMWORD PTR [rcx-32]
	vmaxpd	ymm8, ymm9, ymm3
	vmovapd	ymm3, YMMWORD PTR [rax+64]
	vandpd	ymm4, ymm6, ymm10
	vfmsub132pd	ymm2, ymm0, YMMWORD PTR [rcx+32]
	vmaxpd	ymm6, ymm9, ymm2
	vsubpd	ymm2, ymm1, YMMWORD PTR [rdx+96]
	vfmsub132pd	ymm4, ymm0, YMMWORD PTR [rcx]
	vmovapd	ymm1, YMMWORD PTR [rax+128]
	vmaxpd	ymm7, ymm9, ymm4
	vsubpd	ymm4, ymm3, YMMWORD PTR [rdx+64]
	vandpd	ymm3, ymm2, ymm10
	vsubpd	ymm2, ymm1, YMMWORD PTR [rdx+128]
	vmovapd	ymm1, YMMWORD PTR [rax+160]
	vfmsub132pd	ymm3, ymm0, YMMWORD PTR [rcx+96]
	vandpd	ymm5, ymm4, ymm10
	vmaxpd	ymm4, ymm9, ymm3
	vfmsub132pd	ymm5, ymm0, YMMWORD PTR [rcx+64]
	vmaxpd	ymm5, ymm9, ymm5
	vandpd	ymm3, ymm2, ymm10
	vsubpd	ymm2, ymm1, YMMWORD PTR [rdx+160]
	vfmadd231pd	ymm12, ymm4, ymm4
	vfmsub132pd	ymm3, ymm0, YMMWORD PTR [rcx+128]
	vmaxpd	ymm3, ymm9, ymm3
	vfmadd231pd	ymm11, ymm5, ymm5
	vandpd	ymm1, ymm2, ymm10
	vfmadd231pd	ymm13, ymm3, ymm3
	vfmsub132pd	ymm1, ymm0, YMMWORD PTR [rcx+160]
	vmaxpd	ymm2, ymm9, ymm1
	vmovapd	ymm1, YMMWORD PTR [rax+192]
	vsubpd	ymm1, ymm1, YMMWORD PTR [rdx+192]
	vfmadd231pd	ymm14, ymm2, ymm2
	vandpd	ymm1, ymm1, ymm10
	vfmsub132pd	ymm1, ymm0, YMMWORD PTR [rcx+192]
	vmovapd	ymm0, ymm8
	cmp	rdi, 1024
	vmaxpd	ymm1, ymm9, ymm1
	vfmadd213pd	ymm0, ymm8, YMMWORD PTR [rsp-16]
	vmovapd	ymm8, ymm7
	vmovapd	YMMWORD PTR [rsp-16], ymm0
	vfmadd213pd	ymm8, ymm7, YMMWORD PTR [rsp-48]
	vmovapd	ymm7, ymm6
	vmovapd	YMMWORD PTR [rsp-48], ymm8
	vfmadd213pd	ymm7, ymm6, YMMWORD PTR [rsp-80]
	vfmadd231pd	ymm15, ymm1, ymm1
	vmovapd	YMMWORD PTR [rsp-80], ymm7
	jne	.L6
	vmovapd	ymm6, YMMWORD PTR [rsp-16]
	add	r10, 32
	add	r9, 32
	vmovapd	YMMWORD PTR [r9], ymm8
	vmovapd	YMMWORD PTR [r9-32], ymm6
	vmovapd	YMMWORD PTR [r9+32], ymm7
	vmovapd	YMMWORD PTR [r9+64], ymm11
	vmovapd	YMMWORD PTR [r9+96], ymm12
	vmovapd	YMMWORD PTR [r9+128], ymm13
	vmovapd	YMMWORD PTR [r9+160], ymm14
	vmovapd	YMMWORD PTR [r9+192], ymm15
	cmp	r10, rsi
	jne	.L7
.L5:
	add	r14d, 32
	sub	r15d, 32
	add	r11, 256
	cmp	r14d, DWORD PTR [rsp-84]
	jge	.L11
	mov	edx, DWORD PTR [rsp-84]
	mov	r10, rsi
	sub	edx, r14d
	test	edx, edx
	jg	.L14
.L9:
	mov	rsi, r10
	jmp	.L5
.L11:
	vzeroupper
.L2:
	mov	edx, DWORD PTR [rsp-88]
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
.L8:
	.cfi_restore_state
	xor	esi, esi
	jmp	.L2
	.cfi_endproc
.LFE1081:
	.size	_Z14accumtest4_7_9PdS_S_S_S_Pli, .-_Z14accumtest4_7_9PdS_S_S_S_Pli
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
