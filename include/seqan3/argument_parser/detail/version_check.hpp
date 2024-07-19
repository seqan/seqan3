// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the version check functionality.
 */

#pragma once

#include <array>
#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <optional>
#include <regex>
#include <seqan3/std/charconv>

#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/detail/terminal.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/safe_filesystem_entry.hpp>
#include <seqan3/version.hpp>

#include <sys/stat.h>

namespace seqan3::detail
{

// ---------------------------------------------------------------------------------------------------------------------
// function call_server()
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief Writes a timestamp file and performs the server call to get the newest version information.
 * \param[in] command  The system command as a string. See seqan3::detail::version_checker::command for details.
 * \param[in] prom     A promise object used to track the detached thread which executes this command.
 *
 * This function performs a https server request by executing a hard coded command (string) as a system call.
 */
inline void call_server(std::string const & command, std::promise<bool> prom)
{
    // system call - http response is stored in a file '.config/seqan/{appname}_version'
    if (system(command.c_str()))
        prom.set_value(false);
    else
        prom.set_value(true);
}

// ---------------------------------------------------------------------------------------------------------------------
// version_checker
// ---------------------------------------------------------------------------------------------------------------------

//!\brief A functor whose operator() performs the server http request and version checks.
class version_checker
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief This class has to be initialised with name and version information.
    version_checker() = delete;
    version_checker(version_checker const &) = default;             //!< Defaulted.
    version_checker & operator=(version_checker const &) = default; //!< Defaulted.
    version_checker(version_checker &&) = default;                  //!< Defaulted.
    version_checker & operator=(version_checker &&) = default;      //!< Defaulted.
    ~version_checker() = default;                                   //!< Defaulted.

    /*!\brief Initialises the version_checker with the application name and version.
     * \param[in] name_    The application name.
     * \param[in] version_ The application version.
     * \param[in] app_url  An (github) url with the newest release information of the application.
     */
    version_checker(std::string name_, std::string const & version_, std::string const & app_url = std::string{}) :
        name{std::move(name_)}
    {
        assert(std::regex_match(name, std::regex{"^[a-zA-Z0-9_-]+$"})); // check on construction of the argument parser

        if (!app_url.empty())
        {
            message_app_update.pop_back(); // remove second newline
            message_app_update.append("[APP INFO] :: Visit " + app_url + " for updates.\n\n");
        }

#if defined(NDEBUG)
        timestamp_filename = cookie_path / (name + "_usr.timestamp");
#else
        timestamp_filename = cookie_path / (name + "_dev.timestamp");
#endif
        std::smatch versionMatch;

        // Ensure version string is not corrupt
        if (!version_.empty() && /*regex allows version prefix instead of exact match */
            std::regex_search(version_, versionMatch, std::regex("^([[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).*")))
        {
            version = versionMatch.str(1); // in case the git revision number is given take only version number
        }
    }
    //!\}

    /*!\brief Initialises the version_checker with the application name and version.
     * \param[in] prom The promise to track the state of the detached thread which calls seqan3::detail::call_server.
     *
     * The operator performs the following steps:
     *
     * 1. If a timestamp file already exists, the following behaviour can be expected depending on the file contents:
     *    * content "NEVER": No server call is performed and no update information is printed.
     *    * content "ALWAYS": A server call is performed and update information is printed if a version file exists.
     *    * content is a timestamp: The "time-of-last-version-check" is read from file. The function only
     *                              continues if the last version check is more than a day old.
     *
     * 2. If a version file exists, the app version and seqan3 version are compared to the current ones and the
     *    the following message may be printed:
     *    **Debug mode** (directed at the developer of the application)
     *    * If the app is unregistered (no version information is available at the server) the developer will be
     *      notified that he has the possibility of registering his application with us
     *      (see seqan3::version_checker::message_unregistered_app).
     *    * If the current seqan version is smaller then the one returned by the server call, the developer is notified
     *      that he may update to the newest seqan3 version (see seqan3::version_checker::message_seqan3_update).
     *    * If the current app version is greater than the one returned by the server call, we assume that the
     *      developer has released a new version and is notified to send us the new version
     *      (see seqan3::version_checker::message_registered_app_update).
     *    **Release mode** (directed at the user of the application):
     *    * If the current app version is lower than the one returned by the server call, the user is notified that
     *      a newer version exists.
     */
    void operator()(std::promise<bool> prom)
    {
        std::array<int, 3> empty_version{0, 0, 0};
        std::array<int, 3> srv_app_version{};
        std::array<int, 3> srv_seqan_version{};

        std::ifstream version_file{cookie_path / (name + ".version")};

        if (version_file.is_open())
        {
            std::string line{};
            std::getline(version_file, line); // get first line which should only contain the version number of the app

            if (line != unregistered_app)
                srv_app_version = get_numbers_from_version_string(line);
#if !defined(NDEBUG)
            else
                std::cerr << message_unregistered_app;
#endif // !defined(NDEBUG)

            std::getline(version_file, line); // get second line which should only contain the version number of seqan
            srv_seqan_version = get_numbers_from_version_string(line);

            version_file.close();
        }

#if !defined(NDEBUG) // only check seqan version in debug
        if (srv_seqan_version != empty_version)
        {
            std::array<int, 3> seqan_version = {SEQAN3_VERSION_MAJOR, SEQAN3_VERSION_MINOR, SEQAN3_VERSION_PATCH};

            if (seqan_version < srv_seqan_version)
                std::cerr << message_seqan3_update;
        }
#endif

        if (srv_app_version != empty_version) // app version
        {
#if defined(NDEBUG) // only check app version in release
            if (get_numbers_from_version_string(version) < srv_app_version)
                std::cerr << message_app_update;
#endif // defined(NDEBUG)

#if !defined(NDEBUG) // only notify developer that app version should be updated on server
            if (get_numbers_from_version_string(version) > srv_app_version)
                std::cerr << message_registered_app_update;
#endif // !defined(NDEBUG)
        }

        std::cerr << std::flush;

        std::string program = get_program();

        if (program.empty())
        {
            prom.set_value(false);
            return;
        }

        // 'cookie_path' is no user input and `name` is escaped on construction of the argument parser.
        std::filesystem::path out_file = cookie_path / (name + ".version");

        // build up command for server call
        std::string command = program + // no user defined input
                              " " + out_file.string() + " "
                            + std::string{"http://seqan-update.informatik.uni-tuebingen.de/check/SeqAn3_"} +
#ifdef __linux
                              "Linux" +
#elif __APPLE__
                              "MacOS" +
#elif defined(_WIN32)
                              "Windows" +
#elif __FreeBSD__
                              "FreeBSD" +
#elif __OpenBSD__
                              "OpenBSD" +
#else
                              "unknown" +
#endif
#if __x86_64__ || __ppc64__
                              "_64_" +
#else
                              "_32_" +
#endif
                              name +          // !user input! escaped on construction of the argument parser
                              "_" + version + // !user input! escaped on construction of the version_checker
#if defined(_WIN32)
                              "; exit  [int] -not $?}\" > nul 2>&1";
#else
                              " > /dev/null 2>&1";
#endif

        // launch a separate thread to not defer runtime.
        std::thread(call_server, command, std::move(prom)).detach();
    }

    //!\brief Returns a writable path to store timestamp and version files or an empty path if none exists.
    static std::filesystem::path get_path()
    {
        using namespace std::filesystem;

        path tmp_path;

        tmp_path = std::string{getenv(home_env_name)};
        tmp_path /= ".config";

        // First, create .config if it does not already exist.
        std::error_code err;
        create_directory(tmp_path, err);

        // If this did not fail we, create the seqan subdirectory.
        if (!err)
        {
            tmp_path /= "seqan";
            create_directory(tmp_path, err);
        }

        // .config/seqan cannot be created, try tmp directory.
        if (err)
            tmp_path = temp_directory_path(); // choose temp dir instead

        // check if files can be written inside dir
        path dummy = tmp_path / "dummy.txt";
        std::ofstream file{dummy};
        detail::safe_filesystem_entry file_guard{dummy};

        bool is_open = file.is_open();
        bool is_good = file.good();
        file.close();
        file_guard.remove_no_throw();

        if (!is_good || !is_open) // no write permissions
        {
            tmp_path.clear(); // empty path signals no available directory to write to, version check will not be done
        }

        return tmp_path;
    }

    /*!\brief The central decision whether to perform the version check or not.
     * \param[in] developer_approval Whether the developer approved (update_notifications::on) or not
     *                               (update_notifications::off).
     * \param[in] user_approval      Whether the user approved (true) or not (false) or did not decide (unset optional).
     *
     * The following rules apply:
     *
     * If the developer says no, it rules out all following decisions (even if the user specified --version-check true).
     * No cookie is ever written.
     *
     * If the environment variable SEQAN3_NO_VERSION_CHECK is set no version check is done (rules out all following).
     * No cookie is written.
     *
     * If the user explicitly uses the --version-check option (user_approval is set) it rules out all following
     * decisions. No cookie is written.
     *
     * If none of the above apply, version check was not explicitly handled so cookie content is checked:
     * * NEVER: Do not perform the version check and do not change the cookie.
     * * ALWAYS: Do perform the version check once a day and do not change the cookie.
     * * ASK: Ask the user or default the decision once a day.
     *
     * If the cookie content is "ASK" and the timestamp is older than a day we ask the user,
     * if possible (seqan3::detail::is_terminal()), what he wants to do, set the according cookie for the next time
     * and continue. If we cannot ask the user, the default kicks in (do the check).
     */
    bool decide_if_check_is_performed(update_notifications developer_approval, std::optional<bool> user_approval)
    {
        if (developer_approval == update_notifications::off)
            return false;

        if (std::getenv("SEQAN3_NO_VERSION_CHECK") != nullptr) // environment variable was set
            return false;

        if (user_approval.has_value())
            return user_approval.value();

        // version check was not explicitly handled so let's check the cookie
        if (std::filesystem::exists(cookie_path))
        {
            std::ifstream timestamp_file{timestamp_filename};
            std::string cookie_line{};

            if (timestamp_file.is_open())
            {
                std::getline(timestamp_file, cookie_line); // first line contains the timestamp

                if (get_time_diff_to_current(cookie_line) < 86400 /*one day in seconds*/)
                {
                    return false;
                }

                std::getline(timestamp_file, cookie_line); // second line contains the last user decision

                if (cookie_line == "NEVER")
                {
                    return false;
                }
                else if (cookie_line == "ALWAYS")
                {
                    return true;
                }
                // else we do not return but continue to ask the user

                timestamp_file.close();
            }
        }

        // Up until now, the user did not specify the --version-check option, the environment variable was not set,
        // nor did the the cookie tell us what to do. We will now ask the user if possible or do the check by default.
        write_cookie("ASK"); // Ask again next time when we read the cookie, if this is not overwritten.

        if (detail::is_terminal()) // LCOV_EXCL_START
        {
            std::cerr << R"(
#######################################################################
   Automatic Update Notifications
#######################################################################

 This app can look for updates automatically in the background,
 do you want to do that?

    [a] Always perform version checks for this app (the default).
    [n] Never perform version checks for this app.
    [y] Yes, perform a version check now, and ask again tomorrow.
    [s] Skip the version check now, but ask again tomorrow.

 Please enter one of [a, n, y, s] and press [RETURN].

 For more information, see:
 https://github.com/seqan/seqan3/wiki/Update-Notifications

#######################################################################

)";
            std::string line{};
            std::getline(std::cin, line);
            line.resize(1); // ignore everything but the first char or resizes the empty string to the default

            switch (line[0])
            {
            case 'y':
            {
                return true;
            }
            case 's':
            {
                return false;
            }
            case 'n':
            {
                write_cookie(std::string{"NEVER"}); // overwrite cookie
                return false;
            }
            default:
            {
                write_cookie(std::string{"ALWAYS"}); // overwrite cookie
                return true;
            }
            }
        }
        else // if !detail::is_terminal()
        {
            std::cerr << R"(
#######################################################################
   Automatic Update Notifications
#######################################################################
 This app performs automatic checks for updates. For more information
 see: https://github.com/seqan/seqan3/wiki/Update-Notifications
#######################################################################

)";
            return true; // default: check version if you cannot ask the user
        }
    } // LCOV_EXCL_STOP

    //!\brief The identification string that may appear in the version file if an app is unregistered.
    static constexpr std::string_view unregistered_app = "UNREGISTERED_APP";
    //!\brief The message directed to the developer of the app if a new seqan3 version is available.
    static constexpr std::string_view message_seqan3_update =
        "[SEQAN3 INFO] :: A new SeqAn version is available online.\n"
        "[SEQAN3 INFO] :: Please visit www.github.com/seqan/seqan3.git for an update\n"
        "[SEQAN3 INFO] :: or inform the developer of this app.\n"
        "[SEQAN3 INFO] :: If you don't wish to receive further notifications, set --version-check false.\n\n";
    //!\brief The message directed to the developer of the app if the app is not yet registered with us.
    static constexpr std::string_view message_unregistered_app =
        "[SEQAN3 INFO] :: Thank you for using SeqAn!\n"
        "[SEQAN3 INFO] :: Do you wish to register your app for update notifications?\n"
        "[SEQAN3 INFO] :: Just send an email to support@seqan.de with your app name and version number.\n"
        "[SEQAN3 INFO] :: If you don't wish to receive further notifications, set --version-check false.\n\n";
    //!\brief The message directed to the developer if the application is registered but under a lower version.
    static constexpr std::string_view message_registered_app_update =
        "[APP INFO] :: We noticed the app version you use is newer than the one registered with us.\n"
        "[APP INFO] :: Please send us an email with the new version so we can correct it (support@seqan.de)\n\n";
    //!\brief The message directed to the user of the app if a new app version is available.
    std::string message_app_update =
        "[APP INFO] :: A new version of this application is now available.\n"
        "[APP INFO] :: If you don't wish to receive further notifications, set --version-check false.\n\n";
    /*Might be extended if a url is given on construction.*/

    //!\brief The environment name of the home environment used by getenv()
    static constexpr char const * home_env_name{
#if defined(_WIN32)
        "UserProfile"
#else
        "HOME"
#endif
    };

    //!\brief The application name.
    std::string name;
    //!\brief The version of the application.
    std::string version{"0.0.0"};
    //!\brief The regex to verify a valid version string.
    std::regex version_regex{"^[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+$"};
    //!\brief The path to store timestamp and version files (either ~/.config/seqan or the tmp directory).
    std::filesystem::path cookie_path = get_path();
    //!\brief The timestamp filename.
    std::filesystem::path timestamp_filename;

private:
    //!\brief Returns the command line call as a std::string of an available program depending on the environment.
    static std::string get_program()
    {
#if defined(_WIN32)
        return "powershell.exe -NoLogo -NonInteractive -Command \"& {Invoke-WebRequest -erroraction 'silentlycontinue' "
               "-OutFile";
#else // Unix based platforms.
        if (!system("/usr/bin/env -i wget --version > /dev/null 2>&1"))
            return "/usr/bin/env -i wget --timeout=10 --tries=1 -q -O";
        else if (!system("/usr/bin/env -i curl --version > /dev/null 2>&1"))
            return "/usr/bin/env -i curl --connect-timeout 10 -o";
// In case neither wget nor curl is available try ftp/fetch if system is OpenBSD/FreeBSD.
// Note, both systems have ftp/fetch command installed by default so we do not guard against it.
#    if defined(__OpenBSD__)
        return "/usr/bin/env -i ftp -w10 -Vo";
#    elif defined(__FreeBSD__)
        return "/usr/bin/env -i fetch --timeout=10 -o";
#    else
        return "";
#    endif // __OpenBSD__
#endif     // defined(_WIN32)
    }

    //!\brief Reads the timestamp file if possible and returns the time difference to the current time.
    double get_time_diff_to_current(std::string const & str_time) const
    {
        namespace co = std::chrono;
        double curr = co::duration_cast<co::seconds>(co::system_clock::now().time_since_epoch()).count();

        double d_time{};
        std::from_chars(str_time.data(), str_time.data() + str_time.size(), d_time);

        return curr - d_time;
    }

    /*!\brief Parses a version string into an array of length 3.
     * \param[in] str The version string that must match seqan3::detail::version_regex.
     */
    std::array<int, 3> get_numbers_from_version_string(std::string const & str) const
    {
        std::array<int, 3> result{};

        if (!std::regex_match(str, version_regex))
            return result;

        auto res = std::from_chars(str.data(), str.data() + str.size(), result[0]); // stops and sets res.ptr at '.'
        res = std::from_chars(res.ptr + 1, str.data() + str.size(), result[1]);
        res = std::from_chars(res.ptr + 1, str.data() + str.size(), result[2]);

        return result;
    }

    /*!\brief Writes a cookie file with a specified message.
     * \tparam    msg_type The type of message.
     * \param[in] msg      The message to write into the file (no newline is appended).
     */
    template <typename msg_type>
    void write_cookie(msg_type && msg)
    {
        // The current time
        namespace co = std::chrono;
        auto curr = co::duration_cast<co::seconds>(co::system_clock::now().time_since_epoch()).count();

        std::ofstream timestamp_file{timestamp_filename};

        if (timestamp_file.is_open())
        {
            timestamp_file << curr << '\n' << msg;
            timestamp_file.close();
        }
    }
};

} // namespace seqan3::detail
