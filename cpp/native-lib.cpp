#include <jni.h>
#include <string>
#include <android/log.h>
#include <unistd.h>

#define     TAG                             "bob"

/*
 *              //	Function prototypes
 *
*/

extern "C" long complexQR(double*, double*, double*, int, int, int);
bool throwJavaException(JNIEnv *, std::string, std::string, int, std::string);


extern "C"
JNIEXPORT jlong JNICALL
Java_com_bob_complexqr_MainActivity_complexHouseholder(JNIEnv *env, jobject thiz, jdoubleArray a, jdoubleArray qq, jint rows, jint cols, jint Q) {

    double* aPtr = env->GetDoubleArrayElements(a, nullptr);             // Get C++ pointer to array data
    double* qPtr = env->GetDoubleArrayElements(qq, nullptr);

    double *v;
    if ((v = (double *)calloc(2*rows, sizeof(double))) == nullptr) {
        __android_log_print(ANDROID_LOG_DEBUG, TAG, "malloc failed");
        throwJavaException(env, __FUNCTION__, "malloc failed. OutOfMemoryError", 0, "Exception");
        return 0L;
    }

    __android_log_print(ANDROID_LOG_DEBUG, TAG, "starting complex houseHolder: %d   %d", rows, cols);

//    long time;

    struct timeval  start{}, end{};
    gettimeofday(&start, nullptr);

    long p = complexQR(aPtr, v, qPtr, rows, cols, Q);           // always returns zero

    gettimeofday(&end, nullptr);

    long executionTime = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);

    __android_log_print(ANDROID_LOG_DEBUG, TAG, "houseHolder time, usec: %ld", executionTime);

    env->SetDoubleArrayRegion(a, 0, 2*rows*cols, aPtr);
    env->SetDoubleArrayRegion(qq, 0, 2*rows*rows, qPtr);

    env->ReleaseDoubleArrayElements(a, aPtr, JNI_ABORT);
    env->ReleaseDoubleArrayElements(qq, qPtr, JNI_ABORT);

    free(v);

    return executionTime;

}