
#include <cstdlib>
#include <string>
#include <jni.h>

bool throwJavaException(JNIEnv *env, std::string method_name, std::string exception_msg, int errorCode=0, std::string exception_type= (std::string &) "Exception") {
    char buf[8];
    sprintf(buf, "%d", errorCode);
    std::string code(buf);

    std::string msg = "@" + method_name + ": " + exception_msg + " ";
    if(errorCode != 0) msg += code;

    std::string exception = "java/lang/" + exception_type;

    jclass generalExcep = env->FindClass(exception.c_str());
    if (generalExcep != nullptr) {
        env->ThrowNew(generalExcep, msg.c_str());
        return true;
    }
    return false;
}

bool throwJavaException(JNIEnv *env,  std::string& method_name,  std::string& exception_msg) {
    return throwJavaException(env, method_name, exception_msg, 0, (std::string &)"Exception");
}
