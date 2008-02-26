#ifndef GA_MESSAGE
#define GA_MESSAGE
/* ----------------------------------------------------------------------------
 * R package for graph alignment
 * ----------------------------------------------------------------------------
 *
 * Author: Joern P. Meier <mail@ionflux.org>
 * 
 * The package can be used freely for non-commercial purposes. If you use this 
 * package, the appropriate paper to cite is J. Berg and M. Laessig, 
 * "Cross-species analysis of biological networks by Bayesian alignment", 
 * PNAS 103 (29), 10967-10972 (2006)
 * 
 * This software is made available in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * This software contains code for solving linear assignment problems which was 
 * written by Roy Jonker, MagicLogic Optimization Inc.. Please note that this 
 * code is copyrighted, (c) 2003 MagicLogic Systems Inc., Canada and may be 
 * used for non-commercial purposes only. See 
 * http://www.magiclogic.com/assignment.html for the latest version of the LAP 
 * code and details on licensing.
 *
 * ----------------------------------------------------------------------------
 * Message reporting.
 * ----------------------------------------------------------------------------
 */

/** \file GA_message.h
 * \brief Message reporting.
 *
 * This module provides a simple API for reporting messages. It can be used to 
 * provide information to the user as well as error reporting. The message 
 * handler can be set at run-time, so the right method of output can be 
 * selected easily depending on the project and context in which the module is 
 * being used. A message level can be specified along with the message to 
 * allow a handler to determine the correct method of output. A default 
 * handler is provided for message reporting to the standard output stream. 
 * The default handler also provides some special formatting of messages 
 * depending on the message level.
 */

#ifdef __cplusplus
extern "C"
{
#endif

/** Message level (implementation).
 *
 * The message level is used as additional information for message reporting 
 * functions to determine how the message should be formatted and which method 
 * of output is to be used. Exact behavior of message levels is implementation 
 * defined.
 */
enum GAMessageLevel_Impl
{
    /** Message level: debug.
     */
    GA_MSG_DEBUG = 0,
    /** Message level: info.
     */
    GA_MSG_INFO = 1,
    /** Message level: warning.
     */
    GA_MSG_WARNING = 2,
    /** Message level: error.
     */
    GA_MSG_ERROR = 3
};

/** Message level.
 */
typedef enum GAMessageLevel_Impl GAMessageLevel;

/** Message reporting function.
 *
 * This is the specification for message reporting functions which can then 
 * be set using GA_set_msg_func().
 */
typedef void (*GAMessageFunc)(const char*, GAMessageLevel);

/** Report a message.
 *
 * Report a message. This is the default implementation for message reporting
 * which will write messages to standard output. Some formatting is done 
 * depending on the message level.
 *
 * \param text Text of the message.
 * \param level Message level.
 */
void GA_msg_default(const char* text, GAMessageLevel level);

/** Set message function.
 *
 * Set the message reporting function. This function will then be used for 
 * all subsequent calls to GA_msg().
 *
 * \return Message reporting function.
 *
 * \sa GA_msg
 */
GAMessageFunc GA_set_msg_func(GAMessageFunc msgfunc);

/** Get message function.
 *
 * Get the message reporting function. The message reporting function may be 
 * set with GA_set_msg_func(). This API is a generic, redirectable way of 
 * reporting various informational and debug messages as well as errors. A 
 * default implementation is provided for simple cases. To send a message, use
 *
 * <tt>GA_msg()(&lt;message&gt;, &lt;message level&gt;)</tt>
 *
 * where \c message is the actual message and message level is a value 
 * specifying the type of message. Interpretation of this value is 
 * implementation defined. See GAMessageLevel for values recognized by the 
 * default implementation.
 *
 * \return Message reporting function.
 *
 * \sa GA_set_msg_func
 */
GAMessageFunc GA_msg();

#ifdef __cplusplus
}
#endif
#endif
