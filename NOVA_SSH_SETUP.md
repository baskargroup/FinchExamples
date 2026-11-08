# Let your AI assistant connect to ISU Nova from Ubuntu (so you can ask "status on my Nova runs?")

**For the AI assistant:**
- Read this whole file, then walk the person through it step by step on Ubuntu (or another Linux).
- You make the file edits. The person types every password, verification code and key passphrase,
  in their own terminal. Never ask for, see or store those.
- Replace `<netid>` everywhere with the person's ISU NetID (it is often already in `~/.ssh/config`
  as `User` of a Nova host).

## What this sets up, and why

The normal Nova login asks for a password and a verification code every time. An assistant's tools
have no keyboard, so they can never answer those prompts. ISU's page
https://research.it.iastate.edu/login-nova-ssh-client allows another route:
- SSH key login, but only on `nova-login-1` / `nova-login-2.its.iastate.edu`;
- only from the ISU network (campus network or ISU VPN);
- the key must have a passphrase.

So this setup makes:
- a separate key used only for Nova, protected by a passphrase;
- that key loaded into the key agent of your desktop session, which the assistant's shell can use;
- a shortcut called `nova-key`.

After that, the assistant logs in with no prompt. Your normal login (`ssh <netid>@nova.its.iastate.edu`
with password + code) is unchanged, and any other SSH keys you have keep working as before.

**Needs:**
- Ubuntu with OpenSSH (`ssh -V` works) and `ssh-copy-id` (part of `openssh-client`).
- On the campus network or the ISU VPN.
- `ssh <netid>@nova.its.iastate.edu` already works with password + code.

## Step 0 (assistant): look before changing anything

- `ssh -V`, `which ssh-copy-id ssh-add`.
- `echo $SSH_AUTH_SOCK`: on the Ubuntu desktop this is usually the GNOME keyring agent,
  `/run/user/<uid>/keyring/ssh`. If it is empty, see "No key agent" at the end.
- `cat ~/.ssh/config` (it has no secrets) and check that no `nova-key` entry exists yet.
- Never read private keys; reading `*.pub` files and fingerprints is fine.

## Step 1 (assistant): add the shortcut

1. Copy `~/.ssh/config` (if it exists) to `~/.ssh/config.bak_before_nova_key`.
2. Append these lines. Change nothing else.

```
# Nova, key login for the AI assistant (no prompts)
Host nova-key
    HostName nova-login-1.its.iastate.edu
    User <netid>
    IdentityFile ~/.ssh/id_ed25519_nova
    IdentitiesOnly yes
    BatchMode yes
```

3. Check it: `ssh -G nova-key | grep -E "^(hostname|user|identityfile|batchmode) "`.

Note: if the config has a `Host *` block with its own `IdentityFile`, ssh offers that key too. That is
harmless: Nova rejects it for key login and then accepts the new key.

## Step 2 (the person, in their own terminal)

1. Make the Nova key. Choose a passphrase you will remember; it must not be empty:
   ```bash
   ssh-keygen -t ed25519 -a 100 -f ~/.ssh/id_ed25519_nova -C "<netid> nova-key"
   ```
2. Load it into the key agent. This asks for the passphrase. Repeat it after every reboot or new login:
   ```bash
   ssh-add ~/.ssh/id_ed25519_nova
   ssh-add -l        # should list a line ending in "<netid> nova-key (ED25519)"
   ```
3. Register the public half on Nova with one normal login (password + code):
   ```bash
   ssh-copy-id -i ~/.ssh/id_ed25519_nova.pub <netid>@nova.its.iastate.edu
   ```
   It should end with "Number of key(s) added: 1". It only appends the key to `~/.ssh/authorized_keys`
   on Nova; your home folder is the same on all Nova login nodes.
4. Accept the login node's host key once. Answer `yes` to the fingerprint question; there should be
   no password prompt now:
   ```bash
   ssh -o BatchMode=no nova-key hostname      # prints nova-login-1
   ```

## Step 3 (assistant): test it

```bash
ssh -T -o ConnectTimeout=15 nova-key 'echo NOVA_KEY_OK; hostname; whoami; date'
```

Expected: `NOVA_KEY_OK`, `nova-login-1`, the NetID, the date, with no prompt. The assistant can also
check that it sees the key: `ssh-add -l` lists `<netid> nova-key`.

## "Status on my Nova runs"

Send the commands through a quoted heredoc. The `grep` hides a harmless `stty` line that Nova prints.

```bash
ssh -T nova-key bash -s <<'EOF' 2>&1 | grep -v 'Inappropriate ioctl'
squeue -u $USER -o "%.10i %.20j %.8T %.12M %.12l %.5C %.20R"
sacct -u $USER -S $(date -d '2 days ago' +%F) -X -o JobID,JobName%20,State,Start,End,Elapsed,ExitCode
EOF
```

- Aliases from the person's `~/.bashrc` on Nova (for example `sq`) are not loaded in a plain
  `ssh nova-key cmd`. Use an interactive shell for them: `ssh -T nova-key "bash -ic 'sq'"`.
- To copy files, use `scp` / `rsync` with `nova-key` as the host, for example
  `scp nova-key:/work/<group>/<user>/run/output/log.csv .`
- Keep remote commands simple. Do not nest heredocs inside `bash -s < <(...)`: quoting goes wrong and a
  file can end up overwritten. Copy files with `scp` and run short commands instead.

**Ground rules:**
- Checks like `squeue`, `sacct` and reading logs are fine.
- Ask the person before `sbatch`, `scancel`, uploads or any change on Nova.
- Nova is shared: never write in `/tmp` or in anyone else's folders.
- If the person keeps the code in git: change and commit on the local machine only; on Nova only
  `git pull`.

## Optional: let Nova reach GitHub (for `git clone` / `git pull` on Nova)

Nova has its own SSH keys in `~/.ssh` (for example `id_ed25519`), which GitHub does not know. Without
this, `git clone git@github.com:...` on Nova fails with `Permission denied (publickey)`.

1. On Nova: `cat ~/.ssh/id_ed25519.pub` (if there is none: `ssh-keygen -t ed25519` first, with a
   passphrase).
2. Add that one line in GitHub: Settings > SSH and GPG keys > New SSH key (title e.g. "Nova").
3. Test on Nova: `ssh -T git@github.com` should say "Hi <github-user>! You've successfully authenticated".
4. Clone the branch you need: `git clone -b <branch> git@github.com:<org>/<repo>.git`.

## If it stops working

- Are you on the campus network or the ISU VPN?
- Does `ssh-add -l` list `<netid> nova-key`? After a reboot or new login it does not: run
  `ssh-add ~/.ssh/id_ed25519_nova` again (you type the passphrase).
- If nova-login-1 is down, change the HostName to `nova-login-2.its.iastate.edu` and accept its host
  key once (Step 2.4).
- **No key agent** (`SSH_AUTH_SOCK` empty, e.g. when not logged in to the Ubuntu desktop): start one
  for the session with `eval "$(ssh-agent -s)"`, then `ssh-add ~/.ssh/id_ed25519_nova`. The assistant
  must be started from that same terminal so it inherits `SSH_AUTH_SOCK`. Alternatively add
  `IdentityAgent <path to the agent socket>` to the `nova-key` entry.
- To undo everything: remove the `nova-key` block from `~/.ssh/config` (or restore
  `~/.ssh/config.bak_before_nova_key`), delete `~/.ssh/id_ed25519_nova*`, and remove the line ending in
  `<netid> nova-key` from `~/.ssh/authorized_keys` on Nova.
