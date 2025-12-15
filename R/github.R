library(usethis)
# 1. 告诉 Git 你是谁 (如果你之前没设置过)
use_git_config(user.name = "yantongwan", user.email = "wytsmu@163.com")
# 2. 初始化本地 Git (如果问你是否 commit，选 Yes/Definitely)
use_git()
# 3. 创建 GitHub 令牌 (Token)
# 运行这行后，会自动打开浏览器跳转到 GitHub。
# 在网页上点击 "Generate token"，复制那个以 ghp_ 开头的长字符串。
create_github_token()

edit_r_environ()
# 4. 设置令牌 (运行后，在弹出的窗口里粘贴刚才复制的 Token)
library(gitcreds)
gitcreds_set()
# 5. 一键上传到 GitHub！(见证奇迹的时刻)
usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
use_github()
