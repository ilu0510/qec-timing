; ModuleID = 'hugr-qir'
source_filename = "hugr-qir"
target datalayout = "e-m:e-i8:8:32-i16:16:32-i64:64-i128:128-n32:64-S128"
target triple = "aarch64-unknown-linux-gnu"

%Qubit = type opaque
%Result = type opaque

@0 = private unnamed_addr constant [5 x i8] c"s0_0\00", align 1
@1 = private unnamed_addr constant [5 x i8] c"s0_1\00", align 1
@2 = private unnamed_addr constant [5 x i8] c"s0_2\00", align 1
@3 = private unnamed_addr constant [5 x i8] c"s0_3\00", align 1
@4 = private unnamed_addr constant [5 x i8] c"s1_0\00", align 1
@5 = private unnamed_addr constant [5 x i8] c"s1_1\00", align 1
@6 = private unnamed_addr constant [5 x i8] c"s1_2\00", align 1
@7 = private unnamed_addr constant [5 x i8] c"s1_3\00", align 1
@8 = private unnamed_addr constant [3 x i8] c"m0\00", align 1
@9 = private unnamed_addr constant [3 x i8] c"m1\00", align 1
@10 = private unnamed_addr constant [3 x i8] c"m2\00", align 1
@11 = private unnamed_addr constant [3 x i8] c"m3\00", align 1
@12 = private unnamed_addr constant [3 x i8] c"m4\00", align 1
@13 = private unnamed_addr constant [3 x i8] c"m5\00", align 1
@14 = private unnamed_addr constant [3 x i8] c"m6\00", align 1
@15 = private unnamed_addr constant [3 x i8] c"m7\00", align 1
@16 = private unnamed_addr constant [3 x i8] c"m8\00", align 1
@gen_name = private unnamed_addr constant [8 x i8] c"hugr-qir", section ",qir_generator"
@gen_version = private unnamed_addr constant [5 x i8] c"0.1.2", section ",qir_generator"

define dso_local void @__hugr__.main.1() local_unnamed_addr #0 {
alloca_block:
  tail call void @__quantum__rt__initialize(i8* null)
  %shot = tail call i64 @___get_current_shot()
  %0 = add i64 %shot, 7000000
  tail call void @___random_seed(i64 %0)
  %rintb = tail call i32 @___random_int_bounded(i32 1048576)
  %1 = icmp slt i32 %rintb, 52429
  br i1 %1, label %bb, label %cond_847_case_1

cond_847_case_1:                                  ; preds = %alloca_block, %bb
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* null, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* null)
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb231 = tail call i32 @___random_int_bounded(i32 1048576)
  %2 = icmp slt i32 %rintb231, 52429
  br i1 %2, label %bb0, label %cond_exit_904

bb:                                               ; preds = %alloca_block
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* null)
  br label %cond_847_case_1

bb0:                                              ; preds = %cond_847_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  br label %cond_exit_904

cond_exit_904:                                    ; preds = %bb0, %cond_847_case_1
  %rintb314 = tail call i32 @___random_int_bounded(i32 1048576)
  %3 = icmp slt i32 %rintb314, 52429
  br i1 %3, label %bb1, label %cond_exit_961

bb1:                                              ; preds = %cond_exit_904
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  br label %cond_exit_961

cond_1018_case_1:                                 ; preds = %cond_exit_961, %bb2
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb514 = tail call i32 @___random_int_bounded(i32 1048576)
  %4 = icmp slt i32 %rintb514, 52429
  br i1 %4, label %bb3, label %cond_1384_case_1

cond_exit_961:                                    ; preds = %bb1, %cond_exit_904
  %rintb397 = tail call i32 @___random_int_bounded(i32 1048576)
  %5 = icmp slt i32 %rintb397, 52429
  br i1 %5, label %bb2, label %cond_1018_case_1

bb2:                                              ; preds = %cond_exit_961
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  br label %cond_1018_case_1

cond_1384_case_1:                                 ; preds = %cond_1018_case_1, %bb3
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* null)
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %6 = tail call i1 @__quantum__rt__read_result(%Result* null)
  br i1 %6, label %cond_1366_case_1, label %cond_exit_1366

bb3:                                              ; preds = %cond_1018_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  br label %cond_1384_case_1

cond_1366_case_1:                                 ; preds = %cond_1384_case_1
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_1366

cond_exit_1366:                                   ; preds = %cond_1384_case_1, %cond_1366_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb680 = tail call i32 @___random_int_bounded(i32 1048576)
  %7 = icmp slt i32 %rintb680, 52429
  br i1 %7, label %bb4, label %cond_exit_1132

bb4:                                              ; preds = %cond_exit_1366
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  br label %cond_exit_1132

cond_exit_1132:                                   ; preds = %bb4, %cond_exit_1366
  %rintb763 = tail call i32 @___random_int_bounded(i32 1048576)
  %8 = icmp slt i32 %rintb763, 52429
  br i1 %8, label %bb5, label %cond_exit_1189

bb5:                                              ; preds = %cond_exit_1132
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*))
  br label %cond_exit_1189

cond_1246_case_1:                                 ; preds = %cond_exit_1189, %bb6
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb963 = tail call i32 @___random_int_bounded(i32 1048576)
  %9 = icmp slt i32 %rintb963, 52429
  br i1 %9, label %bb7, label %cond_1483_case_1

cond_exit_1189:                                   ; preds = %bb5, %cond_exit_1132
  %rintb846 = tail call i32 @___random_int_bounded(i32 1048576)
  %10 = icmp slt i32 %rintb846, 52429
  br i1 %10, label %bb6, label %cond_1246_case_1

bb6:                                              ; preds = %cond_exit_1189
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*))
  br label %cond_1246_case_1

cond_1483_case_1:                                 ; preds = %cond_1246_case_1, %bb7
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 1 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %11 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 1 to %Result*))
  br i1 %11, label %cond_1465_case_1, label %cond_exit_1465

bb7:                                              ; preds = %cond_1246_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*))
  br label %cond_1483_case_1

cond_1465_case_1:                                 ; preds = %cond_1483_case_1
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_1465

cond_exit_1465:                                   ; preds = %cond_1483_case_1, %cond_1465_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 2 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %12 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 2 to %Result*))
  br i1 %12, label %cond_1530_case_1, label %cond_exit_1530

cond_1530_case_1:                                 ; preds = %cond_exit_1465
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_1530

cond_exit_1530:                                   ; preds = %cond_exit_1465, %cond_1530_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 3 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %13 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 3 to %Result*))
  br i1 %13, label %cond_1595_case_1, label %cond_exit_1595

cond_1595_case_1:                                 ; preds = %cond_exit_1530
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_1595

cond_exit_1595:                                   ; preds = %cond_exit_1530, %cond_1595_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__rt__bool_record_output(i1 %6, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @0, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %11, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @1, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %12, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @2, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %13, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @3, i64 0, i64 0))
  %rintb1282 = tail call i32 @___random_int_bounded(i32 1048576)
  %14 = icmp slt i32 %rintb1282, 52429
  br i1 %14, label %bb8, label %cond_1688_case_1

cond_1688_case_1:                                 ; preds = %cond_exit_1595, %bb8
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* null, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* null)
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb1382 = tail call i32 @___random_int_bounded(i32 1048576)
  %15 = icmp slt i32 %rintb1382, 52429
  br i1 %15, label %bb9, label %cond_exit_1745

bb8:                                              ; preds = %cond_exit_1595
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* null)
  br label %cond_1688_case_1

bb9:                                              ; preds = %cond_1688_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  br label %cond_exit_1745

cond_exit_1745:                                   ; preds = %bb9, %cond_1688_case_1
  %rintb1465 = tail call i32 @___random_int_bounded(i32 1048576)
  %16 = icmp slt i32 %rintb1465, 52429
  br i1 %16, label %bb10, label %cond_exit_1802

bb10:                                             ; preds = %cond_exit_1745
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  br label %cond_exit_1802

cond_1859_case_1:                                 ; preds = %cond_exit_1802, %bb11
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 1 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb1665 = tail call i32 @___random_int_bounded(i32 1048576)
  %17 = icmp slt i32 %rintb1665, 52429
  br i1 %17, label %bb12, label %cond_2225_case_1

cond_exit_1802:                                   ; preds = %bb10, %cond_exit_1745
  %rintb1548 = tail call i32 @___random_int_bounded(i32 1048576)
  %18 = icmp slt i32 %rintb1548, 52429
  br i1 %18, label %bb11, label %cond_1859_case_1

bb11:                                             ; preds = %cond_exit_1802
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  br label %cond_1859_case_1

cond_2225_case_1:                                 ; preds = %cond_1859_case_1, %bb12
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 4 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %19 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 4 to %Result*))
  br i1 %19, label %cond_2207_case_1, label %cond_exit_2207

bb12:                                             ; preds = %cond_1859_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  br label %cond_2225_case_1

cond_2207_case_1:                                 ; preds = %cond_2225_case_1
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_2207

cond_exit_2207:                                   ; preds = %cond_2225_case_1, %cond_2207_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 4 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb1831 = tail call i32 @___random_int_bounded(i32 1048576)
  %20 = icmp slt i32 %rintb1831, 52429
  br i1 %20, label %bb13, label %cond_exit_1973

bb13:                                             ; preds = %cond_exit_2207
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  br label %cond_exit_1973

cond_exit_1973:                                   ; preds = %bb13, %cond_exit_2207
  %rintb1914 = tail call i32 @___random_int_bounded(i32 1048576)
  %21 = icmp slt i32 %rintb1914, 52429
  br i1 %21, label %bb14, label %cond_exit_2030

bb14:                                             ; preds = %cond_exit_1973
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*))
  br label %cond_exit_2030

cond_2087_case_1:                                 ; preds = %cond_exit_2030, %bb15
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %rintb2114 = tail call i32 @___random_int_bounded(i32 1048576)
  %22 = icmp slt i32 %rintb2114, 52429
  br i1 %22, label %bb16, label %cond_2324_case_1

cond_exit_2030:                                   ; preds = %bb14, %cond_exit_1973
  %rintb1997 = tail call i32 @___random_int_bounded(i32 1048576)
  %23 = icmp slt i32 %rintb1997, 52429
  br i1 %23, label %bb15, label %cond_2087_case_1

bb15:                                             ; preds = %cond_exit_2030
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 7 to %Qubit*))
  br label %cond_2087_case_1

cond_2324_case_1:                                 ; preds = %cond_2087_case_1, %bb16
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 5 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %24 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 5 to %Result*))
  br i1 %24, label %cond_2306_case_1, label %cond_exit_2306

bb16:                                             ; preds = %cond_2087_case_1
  tail call void @__quantum__qis__rz__body(double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 8 to %Qubit*))
  br label %cond_2324_case_1

cond_2306_case_1:                                 ; preds = %cond_2324_case_1
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_2306

cond_exit_2306:                                   ; preds = %cond_2324_case_1, %cond_2306_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 3 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 6 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 6 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %25 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 6 to %Result*))
  br i1 %25, label %cond_2371_case_1, label %cond_exit_2371

cond_2371_case_1:                                 ; preds = %cond_exit_2306
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_2371

cond_exit_2371:                                   ; preds = %cond_exit_2306, %cond_2371_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 2 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0xBFF921FB54442D18, double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rzz__body(double 0x3FF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*), %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 5 to %Qubit*))
  tail call void @__quantum__qis__phasedx__body(double 0x3FF921FB54442D18, double 0x400921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__rz__body(double 0xBFF921FB54442D18, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*), %Result* nonnull inttoptr (i64 7 to %Result*))
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  %26 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 7 to %Result*))
  br i1 %26, label %cond_2436_case_1, label %cond_exit_2436

cond_2436_case_1:                                 ; preds = %cond_exit_2371
  tail call void @__quantum__qis__phasedx__body(double 0x400921FB54442D18, double 0.000000e+00, %Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  br label %cond_exit_2436

cond_exit_2436:                                   ; preds = %cond_exit_2371, %cond_2436_case_1
  tail call void @__quantum__qis__reset__body(%Qubit* nonnull inttoptr (i64 9 to %Qubit*))
  tail call void @__quantum__rt__bool_record_output(i1 %19, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @4, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %24, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @5, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %25, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @6, i64 0, i64 0))
  tail call void @__quantum__rt__bool_record_output(i1 %26, i8* getelementptr inbounds ([5 x i8], [5 x i8]* @7, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* null, %Result* nonnull inttoptr (i64 8 to %Result*))
  %27 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 8 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %27, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @8, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 1 to %Qubit*), %Result* nonnull inttoptr (i64 9 to %Result*))
  %28 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 9 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %28, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @9, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 2 to %Qubit*), %Result* nonnull inttoptr (i64 10 to %Result*))
  %29 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 10 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %29, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @10, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 3 to %Qubit*), %Result* nonnull inttoptr (i64 11 to %Result*))
  %30 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 11 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %30, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @11, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 4 to %Qubit*), %Result* nonnull inttoptr (i64 12 to %Result*))
  %31 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 12 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %31, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @12, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 5 to %Qubit*), %Result* nonnull inttoptr (i64 13 to %Result*))
  %32 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 13 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %32, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @13, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 6 to %Qubit*), %Result* nonnull inttoptr (i64 14 to %Result*))
  %33 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 14 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %33, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @14, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 7 to %Qubit*), %Result* nonnull inttoptr (i64 15 to %Result*))
  %34 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 15 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %34, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @15, i64 0, i64 0))
  tail call void @__quantum__qis__mz__body(%Qubit* nonnull inttoptr (i64 8 to %Qubit*), %Result* nonnull inttoptr (i64 16 to %Result*))
  %35 = tail call i1 @__quantum__rt__read_result(%Result* nonnull inttoptr (i64 16 to %Result*))
  tail call void @__quantum__rt__bool_record_output(i1 %35, i8* getelementptr inbounds ([3 x i8], [3 x i8]* @16, i64 0, i64 0))
  ret void
}

declare i64 @___get_current_shot() local_unnamed_addr

declare void @___random_seed(i64) local_unnamed_addr

declare i32 @___random_int_bounded(i32) local_unnamed_addr

declare void @__quantum__qis__rz__body(double, %Qubit*) local_unnamed_addr

declare void @__quantum__qis__phasedx__body(double, double, %Qubit*) local_unnamed_addr

declare void @__quantum__qis__rzz__body(double, %Qubit*, %Qubit*) local_unnamed_addr

declare void @__quantum__qis__mz__body(%Qubit*, %Result* writeonly) local_unnamed_addr #1

declare void @__quantum__qis__reset__body(%Qubit*) local_unnamed_addr

declare i1 @__quantum__rt__read_result(%Result* readonly) local_unnamed_addr

declare void @__quantum__rt__bool_record_output(i1, i8*) local_unnamed_addr

declare void @__quantum__rt__initialize(i8*) local_unnamed_addr

attributes #0 = { "entry_point" "output_labeling_schema" "qir_profiles"="adaptive_profile" "required_num_qubits"="10" "required_num_results"="17" }
attributes #1 = { "irreversible" }

!llvm.module.flags = !{!0, !1, !2, !3}

!0 = !{i32 1, !"qir_major_version", i32 1}
!1 = !{i32 7, !"qir_minor_version", i32 0}
!2 = !{i32 1, !"dynamic_qubit_management", i1 false}
!3 = !{i32 1, !"dynamic_result_management", i1 false}
