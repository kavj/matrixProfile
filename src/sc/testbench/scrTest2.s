	.file	"scrTest2.c"
	.intel_syntax noprefix
	.text
	.p2align 4,,15
	.globl	sc_init
	.type	sc_init, @function
sc_init:
.LFB10:
	.cfi_startproc
	push	r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	mov	r13d, edx
	push	r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	mov	r12, rdi
	mov	edi, 40
	push	rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	movsx	rbp, esi
	push	rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	sub	rsp, 24
	.cfi_def_cfa_offset 64
	call	malloc
	test	rax, rax
	mov	rbx, rax
	je	.L2
	mov	DWORD PTR [rax+24], ebp
	sal	rbp, 3
	mov	QWORD PTR [rax], r12
	mov	rdi, rbp
	mov	DWORD PTR [rax+28], r13d
	call	malloc
	mov	rdi, rbp
	mov	r12, rax
	mov	QWORD PTR [rbx+8], rax
	call	malloc
	test	r12, r12
	mov	QWORD PTR [rbx+16], rax
	je	.L3
	test	rax, rax
	je	.L3
.L2:
	add	rsp, 24
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	mov	rax, rbx
	pop	rbx
	.cfi_def_cfa_offset 32
	pop	rbp
	.cfi_def_cfa_offset 24
	pop	r12
	.cfi_def_cfa_offset 16
	pop	r13
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L3:
	.cfi_restore_state
	mov	rdi, r12
	mov	QWORD PTR [rsp+8], rax
	call	free
	mov	rax, QWORD PTR [rsp+8]
	mov	rdi, rax
	call	free
	mov	rdi, rbx
	xor	ebx, ebx
	call	free
	jmp	.L2
	.cfi_endproc
.LFE10:
	.size	sc_init, .-sc_init
	.p2align 4,,15
	.globl	mp_init
	.type	mp_init, @function
mp_init:
.LFB11:
	.cfi_startproc
	push	r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	push	r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	push	r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	mov	r12d, esi
	push	rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	mov	ebp, edi
	mov	edi, 16
	push	rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	sub	rsp, 16
	.cfi_def_cfa_offset 64
	call	malloc
	test	rax, rax
	mov	rbx, rax
	je	.L13
	sub	ebp, r12d
	lea	r12d, [rbp+1]
	movsx	r14, r12d
	sal	r14, 3
	mov	rdi, r14
	call	malloc
	mov	rdi, r14
	mov	r13, rax
	mov	QWORD PTR [rbx], rax
	call	malloc
	test	r13, r13
	mov	QWORD PTR [rbx+8], rax
	je	.L14
	test	rax, rax
	je	.L14
.L15:
	test	ebp, ebp
	js	.L13
	mov	rsi, QWORD PTR [rbx+8]
	mov	rcx, QWORD PTR [rbx]
	mov	eax, r12d
	lea	rdx, [rsi+rax*4]
	lea	rax, [rcx+rax*8]
	cmp	rcx, rdx
	setnb	dl
	cmp	rsi, rax
	setnb	al
	or	dl, al
	je	.L16
	cmp	r12d, 24
	jbe	.L16
	mov	r8d, r12d
	pxor	xmm0, xmm0
	shr	r8d, 2
	movsd	xmm2, QWORD PTR .LC1[rip]
	lea	edi, [0+r8*4]
	xor	eax, eax
	xor	edx, edx
.L21:
	add	edx, 1
	movsd	QWORD PTR [rcx+rax*2], xmm2
	movsd	QWORD PTR [rcx+8+rax*2], xmm2
	movsd	QWORD PTR [rcx+16+rax*2], xmm2
	movsd	QWORD PTR [rcx+24+rax*2], xmm2
	movdqu	XMMWORD PTR [rsi+rax], xmm0
	add	rax, 16
	cmp	edx, r8d
	jb	.L21
	cmp	r12d, edi
	je	.L13
	movsd	xmm6, QWORD PTR .LC1[rip]
	movsx	rax, edi
	movsd	QWORD PTR [rcx+rax*8], xmm6
	mov	DWORD PTR [rsi+rax*4], 0
	lea	eax, [rdi+1]
	cmp	ebp, eax
	jl	.L13
	add	edi, 2
	cdqe
	cmp	ebp, edi
	movsd	QWORD PTR [rcx+rax*8], xmm6
	mov	DWORD PTR [rsi+rax*4], 0
	jl	.L13
	movsd	xmm6, QWORD PTR .LC1[rip]
	movsx	rdi, edi
	movsd	QWORD PTR [rcx+rdi*8], xmm6
	mov	DWORD PTR [rsi+rdi*4], 0
.L13:
	add	rsp, 16
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	mov	rax, rbx
	pop	rbx
	.cfi_def_cfa_offset 40
	pop	rbp
	.cfi_def_cfa_offset 32
	pop	r12
	.cfi_def_cfa_offset 24
	pop	r13
	.cfi_def_cfa_offset 16
	pop	r14
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L16:
	.cfi_restore_state
	movsd	xmm1, QWORD PTR .LC1[rip]
	xor	edx, edx
	.p2align 4,,10
	.p2align 3
.L23:
	movsd	QWORD PTR [rcx+rdx*8], xmm1
	mov	DWORD PTR [rsi+rdx*4], 0
	add	rdx, 1
	cmp	ebp, edx
	jge	.L23
	add	rsp, 16
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	mov	rax, rbx
	pop	rbx
	.cfi_def_cfa_offset 40
	pop	rbp
	.cfi_def_cfa_offset 32
	pop	r12
	.cfi_def_cfa_offset 24
	pop	r13
	.cfi_def_cfa_offset 16
	pop	r14
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L14:
	.cfi_restore_state
	mov	rdi, r13
	mov	QWORD PTR [rsp+8], rax
	call	free
	mov	rax, QWORD PTR [rsp+8]
	mov	rdi, rax
	call	free
	mov	rdi, rbx
	xor	ebx, ebx
	call	free
	jmp	.L15
	.cfi_endproc
.LFE11:
	.size	mp_init, .-mp_init
	.p2align 4,,15
	.globl	sc_destroy
	.type	sc_destroy, @function
sc_destroy:
.LFB12:
	.cfi_startproc
	push	rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	mov	rbx, rdi
	mov	rdi, QWORD PTR [rdi+8]
	call	free
	mov	rdi, QWORD PTR [rbx+16]
	call	free
	mov	rdi, rbx
	pop	rbx
	.cfi_def_cfa_offset 8
	jmp	free
	.cfi_endproc
.LFE12:
	.size	sc_destroy, .-sc_destroy
	.p2align 4,,15
	.globl	mp_destroy
	.type	mp_destroy, @function
mp_destroy:
.LFB13:
	.cfi_startproc
	push	rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	mov	rbx, rdi
	mov	rdi, QWORD PTR [rdi]
	call	free
	mov	rdi, QWORD PTR [rbx+8]
	call	free
	mov	rdi, rbx
	pop	rbx
	.cfi_def_cfa_offset 8
	jmp	free
	.cfi_endproc
.LFE13:
	.size	mp_destroy, .-mp_destroy
	.p2align 4,,15
	.globl	iterCnt
	.type	iterCnt, @function
iterCnt:
.LFB14:
	.cfi_startproc
	push	r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	sub	edx, ecx
	mov	r12, rdi
	push	rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	mov	ebp, esi
	push	rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	lea	ebx, [rdx+1]
	add	edx, 2
	imul	edx, ebx
	mov	eax, edx
	shr	eax, 31
	add	edx, eax
	sar	edx
	cvtsi2sd	xmm1, edx
	mulsd	xmm1, xmm0
	movapd	xmm0, xmm1
	call	ceil
	cvttsd2si	ecx, xmm0
	movsx	rdx, ebp
	movsx	rax, ebx
	movsx	rcx, ecx
	add	rdx, rcx
	cmp	rdx, rax
	jle	.L43
	mov	ecx, ebx
	sub	ecx, ebp
	movsx	rcx, ecx
.L43:
	test	ebx, ebx
	jle	.L48
	test	rcx, rcx
	jle	.L48
	sub	ebx, DWORD PTR [r12]
	mov	esi, ebx
	mov	eax, ebx
	jmp	.L47
	.p2align 4,,10
	.p2align 3
.L46:
	add	eax, esi
.L47:
	movsx	rdx, eax
	cmp	rcx, rdx
	jg	.L46
.L49:
	pop	rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	pop	rbp
	.cfi_def_cfa_offset 16
	pop	r12
	.cfi_def_cfa_offset 8
	ret
.L48:
	.cfi_restore_state
	xor	eax, eax
	jmp	.L49
	.cfi_endproc
.LFE14:
	.size	iterCnt, .-iterCnt
	.p2align 4,,15
	.globl	winmeansig
	.type	winmeansig, @function
winmeansig:
.LFB19:
	.cfi_startproc
	push	r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	sub	rsp, 168
	.cfi_def_cfa_offset 224
	cmp	r8d, 1
	mov	QWORD PTR [rsp+8], rdi
	movsd	xmm5, QWORD PTR [rdi]
	jle	.L66
	xorpd	xmm8, xmm8
	lea	rdi, [rdi+8]
	mov	r9d, 1
	movapd	xmm7, xmm8
	jmp	.L53
	.p2align 4,,10
	.p2align 3
.L70:
	mov	r9d, eax
.L53:
	cvtsi2sd	xmm0, r9d
	lea	eax, [r9+1]
	add	rdi, 8
	movsd	xmm1, QWORD PTR [rdi-8]
	cvtsi2sd	xmm2, eax
	cmp	eax, r8d
	subsd	xmm1, xmm5
	mulsd	xmm0, xmm1
	mulsd	xmm0, xmm1
	divsd	xmm1, xmm2
	divsd	xmm0, xmm2
	addsd	xmm5, xmm1
	addsd	xmm7, xmm0
	jne	.L70
.L52:
	cvtsi2sd	xmm4, r8d
	movapd	xmm0, xmm7
	movsd	QWORD PTR [rsi], xmm5
	divsd	xmm0, xmm4
	sqrtsd	xmm1, xmm0
	ucomisd	xmm1, xmm1
	jp	.L71
.L54:
	mov	eax, ecx
	xor	r13d, r13d
	movsx	r14, r8d
	sub	eax, r8d
	mov	r9d, r8d
	mov	DWORD PTR [rsp+24], r8d
	mov	DWORD PTR [rsp], eax
	mov	eax, r8d
	movsd	QWORD PTR [rdx], xmm1
	neg	eax
	lea	r10, [0+r14*8]
	mov	r11, QWORD PTR [rsp+8]
	lea	eax, [rcx+rax*2]
	xor	ecx, ecx
	cmp	r13d, DWORD PTR [rsp]
	lea	rdi, [r14+1]
	mov	DWORD PTR [rsp+16], eax
	lea	eax, [r8-1]
	mov	r8, rcx
	mov	DWORD PTR [rsp+32], eax
	jge	.L72
	.p2align 4,,10
	.p2align 3
.L64:
	cmp	r13d, DWORD PTR [rsp+16]
	jg	.L57
	mov	rax, QWORD PTR [rsp+8]
	cmp	r13d, r9d
	mov	r15d, r9d
	movsd	xmm0, QWORD PTR [rax-8+rdi*8]
	jge	.L73
.L58:
	mov	rcx, rdi
	movapd	xmm1, xmm5
	cvtsi2sd	xmm6, DWORD PTR [rsp+32]
	movapd	xmm3, xmm7
	mov	rbp, r11
	lea	r12d, [r13+1]
	movapd	xmm5, xmm0
	xor	ebx, ebx
	movapd	xmm7, xmm8
	sub	rcx, r8
	jmp	.L65
	.p2align 4,,10
	.p2align 3
.L63:
	cmp	r13d, r12d
	je	.L60
	lea	eax, [rbx+2]
	movsd	xmm2, QWORD PTR [rbp+0+rcx*8]
	cvtsi2sd	xmm9, eax
	lea	eax, [rbx+1]
	subsd	xmm2, xmm5
	cvtsi2sd	xmm0, eax
	mulsd	xmm0, xmm2
	mulsd	xmm0, xmm2
	divsd	xmm2, xmm9
	divsd	xmm0, xmm9
	addsd	xmm5, xmm2
	addsd	xmm7, xmm0
.L60:
	add	rbx, 1
	add	r12d, 1
	add	rbp, 8
.L65:
	movapd	xmm0, xmm1
	movsd	xmm2, QWORD PTR [rbp+0]
	mulsd	xmm0, xmm4
	subsd	xmm0, xmm2
	divsd	xmm0, xmm6
	subsd	xmm2, xmm0
	movapd	xmm1, xmm2
	mulsd	xmm1, xmm6
	mulsd	xmm1, xmm2
	divsd	xmm1, xmm4
	subsd	xmm3, xmm1
	movsd	xmm1, QWORD PTR [rbp+0+r14*8]
	subsd	xmm1, xmm0
	movapd	xmm2, xmm1
	mulsd	xmm2, xmm6
	mulsd	xmm2, xmm1
	divsd	xmm1, xmm4
	divsd	xmm2, xmm4
	addsd	xmm1, xmm0
	movsd	QWORD PTR [rsi+8+rbx*8], xmm1
	addsd	xmm3, xmm2
	movapd	xmm0, xmm3
	divsd	xmm0, xmm4
	sqrtsd	xmm2, xmm0
	ucomisd	xmm2, xmm2
	jp	.L74
.L61:
	cmp	r15d, r12d
	movsd	QWORD PTR [rdx+8+rbx*8], xmm2
	jg	.L63
.L59:
	mov	eax, DWORD PTR [rsp+24]
	add	r11, r10
	add	rsi, r10
	add	rdx, r10
	add	rdi, r14
	add	r8, r14
	add	r13d, eax
	add	r9d, eax
	cmp	r13d, DWORD PTR [rsp]
	jl	.L64
.L72:
	add	rsp, 168
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	pop	rbx
	.cfi_def_cfa_offset 48
	pop	rbp
	.cfi_def_cfa_offset 40
	pop	r12
	.cfi_def_cfa_offset 32
	pop	r13
	.cfi_def_cfa_offset 24
	pop	r14
	.cfi_def_cfa_offset 16
	pop	r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L57:
	.cfi_restore_state
	mov	rax, QWORD PTR [rsp+8]
	mov	r15d, DWORD PTR [rsp]
	movsd	xmm0, QWORD PTR [rax-8+rdi*8]
	jmp	.L58
.L73:
	movapd	xmm5, xmm0
	movapd	xmm7, xmm8
	jmp	.L59
.L66:
	xorpd	xmm8, xmm8
	movapd	xmm7, xmm8
	jmp	.L52
.L74:
	mov	QWORD PTR [rsp+152], rcx
	mov	QWORD PTR [rsp+144], rdx
	movsd	QWORD PTR [rsp+88], xmm8
	mov	QWORD PTR [rsp+136], rdi
	mov	QWORD PTR [rsp+128], r8
	movsd	QWORD PTR [rsp+80], xmm6
	mov	QWORD PTR [rsp+120], rsi
	movsd	QWORD PTR [rsp+72], xmm1
	mov	QWORD PTR [rsp+112], r10
	movsd	QWORD PTR [rsp+64], xmm3
	mov	QWORD PTR [rsp+104], r11
	movsd	QWORD PTR [rsp+56], xmm4
	mov	DWORD PTR [rsp+100], r9d
	movsd	QWORD PTR [rsp+48], xmm5
	movsd	QWORD PTR [rsp+40], xmm7
	call	sqrt
	mov	rcx, QWORD PTR [rsp+152]
	movapd	xmm2, xmm0
	mov	rdx, QWORD PTR [rsp+144]
	mov	rdi, QWORD PTR [rsp+136]
	mov	r8, QWORD PTR [rsp+128]
	mov	rsi, QWORD PTR [rsp+120]
	mov	r10, QWORD PTR [rsp+112]
	mov	r11, QWORD PTR [rsp+104]
	mov	r9d, DWORD PTR [rsp+100]
	movsd	xmm8, QWORD PTR [rsp+88]
	movsd	xmm6, QWORD PTR [rsp+80]
	movsd	xmm1, QWORD PTR [rsp+72]
	movsd	xmm3, QWORD PTR [rsp+64]
	movsd	xmm4, QWORD PTR [rsp+56]
	movsd	xmm5, QWORD PTR [rsp+48]
	movsd	xmm7, QWORD PTR [rsp+40]
	jmp	.L61
.L71:
	movsd	QWORD PTR [rsp+64], xmm8
	mov	DWORD PTR [rsp+56], r8d
	mov	DWORD PTR [rsp+48], ecx
	mov	QWORD PTR [rsp+40], rdx
	mov	QWORD PTR [rsp+32], rsi
	movsd	QWORD PTR [rsp+24], xmm4
	movsd	QWORD PTR [rsp+16], xmm5
	movsd	QWORD PTR [rsp], xmm7
	call	sqrt
	movsd	xmm8, QWORD PTR [rsp+64]
	movapd	xmm1, xmm0
	mov	r8d, DWORD PTR [rsp+56]
	mov	ecx, DWORD PTR [rsp+48]
	mov	rdx, QWORD PTR [rsp+40]
	mov	rsi, QWORD PTR [rsp+32]
	movsd	xmm4, QWORD PTR [rsp+24]
	movsd	xmm5, QWORD PTR [rsp+16]
	movsd	xmm7, QWORD PTR [rsp]
	jmp	.L54
	.cfi_endproc
.LFE19:
	.size	winmeansig, .-winmeansig
	.p2align 4,,15
	.globl	sccomp
	.type	sccomp, @function
sccomp:
.LFB20:
	.cfi_startproc
	push	r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	sub	rsp, 152
	.cfi_def_cfa_offset 208
	mov	ebp, DWORD PTR [rsp+216]
	mov	ebx, DWORD PTR [rsp+224]
	mov	edx, DWORD PTR [rsp+208]
	lea	r15d, [rbp+0+rbx]
	movsx	r12, ebp
	add	edx, ebp
	lea	r13, [0+r12*8]
	movsx	r14, r15d
	cmp	ebp, edx
	movsd	xmm5, QWORD PTR [rsi+r12*8]
	lea	r11, [0+r14*8]
	movsd	xmm6, QWORD PTR [rsi+r14*8]
	jge	.L94
	xorpd	xmm2, xmm2
	sub	edx, ebp
	lea	rax, [rdi+r13]
	sub	edx, 1
	add	rdx, r12
	movapd	xmm4, xmm2
	movapd	xmm3, xmm2
	lea	r10, [rdi+8+rdx*8]
	movsx	rdx, ebx
	.p2align 4,,10
	.p2align 3
.L78:
	movsd	xmm1, QWORD PTR [rax]
	movsd	xmm0, QWORD PTR [rax+rdx*8]
	add	rax, 8
	subsd	xmm1, xmm5
	cmp	rax, r10
	subsd	xmm0, xmm6
	movapd	xmm7, xmm1
	mulsd	xmm7, xmm1
	addsd	xmm3, xmm7
	movapd	xmm7, xmm0
	mulsd	xmm7, xmm0
	mulsd	xmm0, xmm1
	addsd	xmm4, xmm7
	addsd	xmm2, xmm0
	jne	.L78
	movapd	xmm0, xmm3
	mulsd	xmm0, xmm4
.L76:
	sqrtsd	xmm1, xmm0
	ucomisd	xmm1, xmm1
	jp	.L100
.L79:
	movapd	xmm0, xmm2
	lea	rax, [rcx+r13]
	divsd	xmm0, xmm1
	movsd	xmm1, QWORD PTR [rax]
	ucomisd	xmm1, xmm0
	jbe	.L81
	movsd	QWORD PTR [rax], xmm0
	mov	DWORD PTR [r8+r12*4], r15d
.L81:
	add	r11, rcx
	movsd	xmm1, QWORD PTR [r11]
	ucomisd	xmm1, xmm0
	jbe	.L83
	movsd	QWORD PTR [r11], xmm0
	mov	DWORD PTR [r8+r14*4], ebp
.L83:
	sub	r9d, DWORD PTR [rsp+208]
	lea	r15d, [rbp+1]
	cmp	r15d, r9d
	jge	.L75
	movsx	rdx, r15d
	movsx	r11, DWORD PTR [rsp+208]
	lea	rax, [rsi+8+r13]
	lea	r14, [rcx+rdx*8]
	movsx	rcx, ebx
	add	r13, rdi
	mov	rdi, rdx
	lea	rdx, [rcx+rdx]
	add	ebp, ebp
	mov	r10, rcx
	movsx	rbp, ebp
	mov	rsi, r12
	lea	rdx, [r8+rdx*4]
	mov	r8, rdi
	add	r12, r11
	neg	r8
	neg	r10
	add	r11, rbp
	sal	r8, 3
	sal	rsi, 4
	sal	r10, 2
	sub	r12, rdi
	sub	r11, rdi
	jmp	.L93
	.p2align 4,,10
	.p2align 3
.L102:
	movapd	xmm6, xmm8
	movapd	xmm5, xmm7
.L93:
	movsd	xmm11, QWORD PTR [r13+0]
	lea	rdi, [r13+0+r8]
	lea	ebp, [r15+rbx]
	movsd	xmm1, QWORD PTR [r13+8+r12*8]
	movapd	xmm9, xmm11
	movsd	xmm7, QWORD PTR [rax]
	movapd	xmm10, xmm1
	subsd	xmm9, xmm5
	movsd	xmm12, QWORD PTR [rdi+8+rsi]
	subsd	xmm10, xmm5
	movapd	xmm5, xmm11
	subsd	xmm1, xmm7
	movsd	xmm0, QWORD PTR [r13+8+r11*8]
	subsd	xmm5, xmm7
	movsd	xmm8, QWORD PTR [rax+rcx*8]
	movapd	xmm11, xmm0
	subsd	xmm0, xmm6
	mulsd	xmm1, xmm10
	subsd	xmm11, xmm8
	mulsd	xmm5, xmm9
	mulsd	xmm0, xmm11
	subsd	xmm1, xmm5
	movapd	xmm5, xmm12
	subsd	xmm5, xmm8
	addsd	xmm3, xmm1
	movapd	xmm1, xmm12
	subsd	xmm1, xmm6
	mulsd	xmm1, xmm5
	subsd	xmm0, xmm1
	addsd	xmm4, xmm0
	movapd	xmm0, xmm3
	mulsd	xmm0, xmm4
	sqrtsd	xmm1, xmm0
	ucomisd	xmm1, xmm1
	jp	.L101
.L86:
	movapd	xmm0, xmm10
	mulsd	xmm5, xmm9
	mulsd	xmm0, xmm11
	subsd	xmm0, xmm5
	addsd	xmm2, xmm0
	movapd	xmm0, xmm2
	divsd	xmm0, xmm1
	movsd	xmm1, QWORD PTR [r14]
	ucomisd	xmm1, xmm0
	jbe	.L88
	movsd	QWORD PTR [r14], xmm0
	mov	DWORD PTR [rdx+r10], ebp
.L88:
	movsd	xmm1, QWORD PTR [r14+rcx*8]
	ucomisd	xmm1, xmm0
	jbe	.L90
	movsd	QWORD PTR [r14+rcx*8], xmm0
	mov	DWORD PTR [rdx], r15d
.L90:
	add	r15d, 1
	add	rax, 8
	add	r13, 8
	add	r14, 8
	add	rdx, 4
	cmp	r15d, r9d
	jne	.L102
.L75:
	add	rsp, 152
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	pop	rbx
	.cfi_def_cfa_offset 48
	pop	rbp
	.cfi_def_cfa_offset 40
	pop	r12
	.cfi_def_cfa_offset 32
	pop	r13
	.cfi_def_cfa_offset 24
	pop	r14
	.cfi_def_cfa_offset 16
	pop	r15
	.cfi_def_cfa_offset 8
	ret
.L94:
	.cfi_restore_state
	xorpd	xmm0, xmm0
	movapd	xmm2, xmm0
	movapd	xmm4, xmm0
	movapd	xmm3, xmm0
	jmp	.L76
.L101:
	mov	QWORD PTR [rsp+136], r11
	mov	QWORD PTR [rsp+128], r10
	movsd	QWORD PTR [rsp+72], xmm5
	mov	QWORD PTR [rsp+120], r8
	mov	QWORD PTR [rsp+112], rdx
	movsd	QWORD PTR [rsp+64], xmm11
	mov	QWORD PTR [rsp+104], rcx
	movsd	QWORD PTR [rsp+56], xmm9
	mov	QWORD PTR [rsp+96], rsi
	movsd	QWORD PTR [rsp+48], xmm10
	mov	DWORD PTR [rsp+92], r9d
	movsd	QWORD PTR [rsp+40], xmm7
	mov	QWORD PTR [rsp+80], rax
	movsd	QWORD PTR [rsp+32], xmm8
	movsd	QWORD PTR [rsp+24], xmm2
	movsd	QWORD PTR [rsp+16], xmm4
	movsd	QWORD PTR [rsp+8], xmm3
	call	sqrt
	mov	r11, QWORD PTR [rsp+136]
	movapd	xmm1, xmm0
	mov	r10, QWORD PTR [rsp+128]
	mov	r8, QWORD PTR [rsp+120]
	mov	rdx, QWORD PTR [rsp+112]
	mov	rcx, QWORD PTR [rsp+104]
	mov	rsi, QWORD PTR [rsp+96]
	mov	r9d, DWORD PTR [rsp+92]
	mov	rax, QWORD PTR [rsp+80]
	movsd	xmm5, QWORD PTR [rsp+72]
	movsd	xmm11, QWORD PTR [rsp+64]
	movsd	xmm9, QWORD PTR [rsp+56]
	movsd	xmm10, QWORD PTR [rsp+48]
	movsd	xmm7, QWORD PTR [rsp+40]
	movsd	xmm8, QWORD PTR [rsp+32]
	movsd	xmm2, QWORD PTR [rsp+24]
	movsd	xmm4, QWORD PTR [rsp+16]
	movsd	xmm3, QWORD PTR [rsp+8]
	jmp	.L86
.L100:
	mov	DWORD PTR [rsp+92], r9d
	mov	QWORD PTR [rsp+80], r8
	movsd	QWORD PTR [rsp+48], xmm2
	mov	QWORD PTR [rsp+72], rcx
	mov	QWORD PTR [rsp+64], rsi
	movsd	QWORD PTR [rsp+40], xmm4
	mov	QWORD PTR [rsp+56], rdi
	movsd	QWORD PTR [rsp+32], xmm3
	movsd	QWORD PTR [rsp+24], xmm6
	mov	QWORD PTR [rsp+16], r11
	movsd	QWORD PTR [rsp+8], xmm5
	call	sqrt
	mov	r9d, DWORD PTR [rsp+92]
	movapd	xmm1, xmm0
	mov	r8, QWORD PTR [rsp+80]
	mov	rcx, QWORD PTR [rsp+72]
	mov	rsi, QWORD PTR [rsp+64]
	mov	rdi, QWORD PTR [rsp+56]
	movsd	xmm2, QWORD PTR [rsp+48]
	mov	r11, QWORD PTR [rsp+16]
	movsd	xmm4, QWORD PTR [rsp+40]
	movsd	xmm3, QWORD PTR [rsp+32]
	movsd	xmm6, QWORD PTR [rsp+24]
	movsd	xmm5, QWORD PTR [rsp+8]
	jmp	.L79
	.cfi_endproc
.LFE20:
	.size	sccomp, .-sccomp
	.p2align 4,,15
	.globl	corrToDist
	.type	corrToDist, @function
corrToDist:
.LFB22:
	.cfi_startproc
	sub	esi, edx
	js	.L111
	add	edx, edx
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	mov	ebp, esi
	cvtsi2sd	xmm3, edx
	push	rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	xor	ebx, ebx
	sub	rsp, 40
	.cfi_def_cfa_offset 64
	movsd	xmm2, QWORD PTR .LC3[rip]
	.p2align 4,,10
	.p2align 3
.L108:
	movapd	xmm0, xmm2
	subsd	xmm0, QWORD PTR [rdi+rbx*8]
	mulsd	xmm0, xmm3
	sqrtsd	xmm1, xmm0
	ucomisd	xmm1, xmm1
	jp	.L112
.L105:
	movsd	QWORD PTR [rdi+rbx*8], xmm1
	add	rbx, 1
	cmp	ebp, ebx
	jge	.L108
	add	rsp, 40
	.cfi_def_cfa_offset 24
	pop	rbx
	.cfi_restore 3
	.cfi_def_cfa_offset 16
	pop	rbp
	.cfi_restore 6
	.cfi_def_cfa_offset 8
.L111:
	rep ret
.L112:
	.cfi_def_cfa_offset 64
	.cfi_offset 3, -24
	.cfi_offset 6, -16
	movsd	QWORD PTR [rsp+24], xmm2
	mov	QWORD PTR [rsp+16], rdi
	movsd	QWORD PTR [rsp+8], xmm3
	call	sqrt
	movsd	xmm2, QWORD PTR [rsp+24]
	movapd	xmm1, xmm0
	mov	rdi, QWORD PTR [rsp+16]
	movsd	xmm3, QWORD PTR [rsp+8]
	jmp	.L105
	.cfi_endproc
.LFE22:
	.size	corrToDist, .-corrToDist
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC4:
	.string	"base n: %d\n"
	.text
	.p2align 4,,15
	.globl	scBlockSolver
	.type	scBlockSolver, @function
scBlockSolver:
.LFB24:
	.cfi_startproc
	push	r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	mov	r15, rdi
	push	r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	mov	rbx, rsi
	sub	rsp, 40
	.cfi_def_cfa_offset 96
	mov	eax, DWORD PTR [rdi+24]
	mov	r13d, DWORD PTR [rdi+28]
	mov	edi, OFFSET FLAT:.LC4
	mov	r12d, eax
	mov	DWORD PTR [rsp+28], eax
	mov	esi, eax
	xor	eax, eax
	call	printf
	test	r12d, r12d
	jle	.L113
	sub	r12d, r13d
	xor	ebp, ebp
	.p2align 4,,10
	.p2align 3
.L115:
	mov	r14d, DWORD PTR [rsp+28]
	sub	r14d, ebp
	cmp	r13d, r12d
	jge	.L117
	cmp	r14d, 4096
	mov	eax, 4096
	cmovg	r14d, eax
	mov	DWORD PTR [rsp+24], r14d
	mov	r14d, r13d
	.p2align 4,,10
	.p2align 3
.L118:
	mov	r8, QWORD PTR [rbx+8]
	mov	rdx, QWORD PTR [r15+16]
	mov	rsi, QWORD PTR [r15+8]
	mov	r9d, DWORD PTR [rsp+24]
	mov	rcx, QWORD PTR [rbx]
	mov	rdi, QWORD PTR [r15]
	mov	DWORD PTR [rsp+16], r14d
	mov	DWORD PTR [rsp+8], ebp
	add	r14d, 1
	mov	DWORD PTR [rsp], r13d
	call	sccomp
	cmp	r14d, r12d
	jne	.L118
.L117:
	add	ebp, 4096
	sub	r12d, 4096
	cmp	DWORD PTR [rsp+28], ebp
	jg	.L115
.L113:
	add	rsp, 40
	.cfi_def_cfa_offset 56
	pop	rbx
	.cfi_def_cfa_offset 48
	pop	rbp
	.cfi_def_cfa_offset 40
	pop	r12
	.cfi_def_cfa_offset 32
	pop	r13
	.cfi_def_cfa_offset 24
	pop	r14
	.cfi_def_cfa_offset 16
	pop	r15
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE24:
	.size	scBlockSolver, .-scBlockSolver
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	0
	.long	-1074790400
	.align 8
.LC3:
	.long	0
	.long	1072693248
	.ident	"GCC: (GNU) 4.8.5 20150623 (Red Hat 4.8.5-11)"
	.section	.note.GNU-stack,"",@progbits
