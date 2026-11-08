# Let your AI assistant connect to ISU Nova from a Mac (so you can ask "status on my Nova runs?")

macOS version of `NOVA_SSH_SETUP.md`. Same idea, three Mac differences:
1. macOS already runs an SSH key agent for your login session (`SSH_AUTH_SOCK` points at
   `/var/run/com.apple.launchd.*/Listeners`). Nothing to start.
2. The key passphrase can be stored in the macOS Keychain, so the key loads itself after a reboot.
   You never have to run `ssh-add` again.
3. Off campus you must be on the ISU VPN first (Cisco Secure Client). The assistant cannot connect
   the VPN for you: that needs your password and a multi-factor approval.

**For the AI assistant:**
- Read this whole file, then walk the person through it step by step.
- You make the file edits. The person types every password, verification code and key passphrase,
  in their own terminal. Never ask for, see or store those.
- Replace `<netid>` everywhere with the person's ISU NetID.

## What this sets up, and why

The normal Nova login asks for a verification code (TOTP from an authenticator app) and a password
every time. An assistant's tools have no keyboard, so they can never answer those prompts. ISU's page
https://research.it.iastate.edu/login-nova-ssh-client allows another route:
- SSH key login, but only on `nova-login-1` / `nova-login-2.its.iastate.edu`
  (`nova.its.iastate.edu` is a pool, not one server, so key login does not work there);
- only from the ISU network: campus Ethernet, Eduroam Wi-Fi, or the ISU VPN. ISU's Nova page says
  the VPN (or Eduroam) is needed for key login "regardless of whether you are on or off campus";
- the key must have a passphrase;
- TOTP must already be set up, because registering the key needs one normal login.

So this setup makes:
- a separate key used only for Nova, protected by a passphrase, with the passphrase in the Keychain;
- a shortcut called `nova-key`.

After that, the assistant logs in with no prompt. Your normal login (`ssh <netid>@nova.its.iastate.edu`
with code + password) is unchanged, and other SSH keys (for example your GitHub key) keep working.

**Needs:**
- macOS with OpenSSH (`ssh -V`; macOS 26 ships OpenSSH 10) and `ssh-copy-id` (`/usr/bin/ssh-copy-id`).
- On campus / Eduroam, or connected to the ISU VPN (see "Off campus: the ISU VPN" below).
- `ssh <netid>@nova.its.iastate.edu` already works with code + password.

## Off campus: the ISU VPN

Off campus, Nova is simply unreachable: `nova-login-1.its.iastate.edu` does not even resolve in
public DNS, and port 22 on `nova.its.iastate.edu` is closed. The VPN puts the Mac on the ISU network.

One-time install (the person does this):
1. In a browser go to https://vpn.iastate.edu and sign in with your NetID and password
   (plus the multi-factor prompt you use for login.iastate.edu). Download the **Cisco Secure Client**
   for macOS from there, or follow ISU's article "Virtual Private Network: Install and Connect to VPN"
   (KB0011105, needs ISU login):
   https://iastate.service-now.com/it?id=kb_article&sysparm_article=KB0011105
2. Open the `.dmg`, run the installer. macOS asks to allow the Cisco system extension in
   System Settings > General > Login Items & Extensions (or Privacy & Security). Allow it.

Every time you work off campus:
1. Open **Cisco Secure Client** (in /Applications), server `vpn.iastate.edu`, click Connect.
2. Sign in with NetID + password and approve the multi-factor prompt.
3. Only then run `ssh ... nova-key` or ask the assistant for Nova status.

Checks the assistant can run (no credentials involved):
```bash
scutil --nc list                                   # Cisco VPN shows here when installed
ifconfig | grep -A2 utun | grep "inet "            # a utun interface with an address = VPN tunnel up
nc -z -w 5 nova-login-1.its.iastate.edu 22 && echo reachable
```
If the last line is not `reachable`, the VPN is not connected. Ask the person to connect it; do not
retry in a loop.

Alternatives without the VPN: Nova OnDemand in the browser (https://nova-ondemand.its.iastate.edu,
"Clusters > Nova Shell") for you by hand; nothing for the assistant.

## Step 0 (assistant): look before changing anything

- `ssh -V`, `which ssh-copy-id ssh-add`.
- `echo $SSH_AUTH_SOCK` should be `/var/run/com.apple.launchd.<id>/Listeners`. `ssh-add -l` lists
  what is loaded ("The agent has no identities." is fine at this point).
- `cat ~/.ssh/config` (it has no secrets) and check that no `nova-key` entry exists yet.
- Never read private keys; reading `*.pub` files and fingerprints is fine.

## Step 1 (assistant): add the shortcut

1. Copy `~/.ssh/config` (if it exists) to `~/.ssh/config.bak_before_nova_key`.
2. Insert these lines **above any `Host *` block** (ssh takes the first matching value, so the
   `nova-key` block must come first). Change nothing else.

```
# Nova, key login for the AI assistant (no prompts)
Host nova-key
    HostName nova-login-1.its.iastate.edu
    User <netid>
    IdentityFile ~/.ssh/id_ed25519_nova
    IdentitiesOnly yes
    AddKeysToAgent yes
    UseKeychain yes
    BatchMode yes
```

`AddKeysToAgent yes` + `UseKeychain yes` is the macOS part: the first time the key is used, ssh
fetches the passphrase from the Keychain and loads the key into the agent, with no prompt.

3. Check it: `ssh -G nova-key | grep -E "^(hostname|user|identityfile|batchmode|addkeystoagent) "`.

Notes:
- If the config has a `Host *` block with `ControlMaster auto` / `ControlPersist`, that is fine and
  even helpful: the first login stays open for a while and later commands reuse it.
- If the `Host *` block has its own `IdentityFile`, ssh offers that key too. Harmless: Nova rejects
  it and then accepts the new key.

## Step 2 (the person, in their own terminal, on campus / Eduroam / VPN)

1. Make the Nova key. Choose a passphrase you will remember; it must not be empty:
   ```bash
   ssh-keygen -t ed25519 -a 100 -f ~/.ssh/id_ed25519_nova -C "<netid> nova-key"
   ```
2. Load it into the agent and store the passphrase in the Keychain. Asks for the passphrase once:
   ```bash
   ssh-add --apple-use-keychain ~/.ssh/id_ed25519_nova
   ssh-add -l        # should list a line ending in "<netid> nova-key (ED25519)"
   ```
3. Register the public half on Nova with one normal login (verification code, then password):
   ```bash
   ssh-copy-id -i ~/.ssh/id_ed25519_nova.pub <netid>@nova.its.iastate.edu
   ```
   It should end with "Number of key(s) added: 1". It only appends the key to `~/.ssh/authorized_keys`
   on Nova; your home folder is the same on all Nova login nodes.
4. Accept the login node's host key once. Answer `yes` to the fingerprint question (compare with
   https://research.it.iastate.edu/nova-ssh-host-key-fingerprints if you like); there should be no
   password prompt now:
   ```bash
   ssh -o BatchMode=no nova-key hostname      # prints nova-login-1
   ```

## Step 3 (assistant): test it

```bash
ssh -T -o ConnectTimeout=15 nova-key 'echo NOVA_KEY_OK; hostname; whoami; date'
```

Expected: `NOVA_KEY_OK`, `nova-login-1`, the NetID, the date, with no prompt. `ssh-add -l` should list
`<netid> nova-key`. After a reboot the agent is empty again, but the first `ssh nova-key` reloads the
key from the Keychain by itself (that is what `AddKeysToAgent` + `UseKeychain` do).

## "Status on my Nova runs"

Send the commands through a quoted heredoc. The `grep` hides a harmless `stty` line that Nova prints.
The commands run on Nova (Linux), so GNU `date -d` is correct there even though the Mac's `date`
does not have it.

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
- If the person keeps the code in git: change and commit on the Mac only; on Nova only `git pull`.

## Optional: let Nova reach GitHub (for `git clone` / `git pull` on Nova)

Nova has its own SSH keys in `~/.ssh` (for example `id_ed25519`), which GitHub does not know. Without
this, `git clone git@github.com:...` on Nova fails with `Permission denied (publickey)`.

1. On Nova: `cat ~/.ssh/id_ed25519.pub` (if there is none: `ssh-keygen -t ed25519` first, with a
   passphrase).
2. Add that one line in GitHub: Settings > SSH and GPG keys > New SSH key (title e.g. "Nova").
3. Test on Nova: `ssh -T git@github.com` should say "Hi <github-user>! You've successfully authenticated".
4. Clone the branch you need: `git clone -b <branch> git@github.com:<org>/<repo>.git`.

## If it stops working

- Off campus: is the ISU VPN connected? `nc -z -w 5 nova-login-1.its.iastate.edu 22` must succeed.
  On campus Wi-Fi other than Eduroam (for example a guest network) the VPN is needed too.
- Does `ssh-add -l` list `<netid> nova-key`? If not, `ssh-add --apple-use-keychain ~/.ssh/id_ed25519_nova`
  (no passphrase asked when it is in the Keychain; otherwise type it once more).
- If nova-login-1 is down, change the HostName to `nova-login-2.its.iastate.edu` and accept its host
  key once (Step 2.4).
- Assistant started from a terminal with an empty `SSH_AUTH_SOCK` (rare on a Mac, e.g. inside an
  ssh session into the Mac): `eval "$(ssh-agent -s)"`, then `ssh-add --apple-use-keychain ~/.ssh/id_ed25519_nova`,
  and start the assistant from that same terminal.
- A stale ControlMaster socket (`/tmp/ssh-<netid>@nova-login-1.its.iastate.edu:22`) after a VPN drop:
  `ssh -O exit nova-key`, then connect again.
- To undo everything: remove the `nova-key` block from `~/.ssh/config` (or restore
  `~/.ssh/config.bak_before_nova_key`), delete `~/.ssh/id_ed25519_nova*`, remove the passphrase
  from Keychain Access (search "SSH: id_ed25519_nova"), and remove the line ending in
  `<netid> nova-key` from `~/.ssh/authorized_keys` on Nova.

## Sources
- https://research.it.iastate.edu/login-nova-ssh-client (key login, login nodes, VPN/Eduroam rule)
- https://research.it.iastate.edu/nova (MFA/TOTP setup, OnDemand)
- https://research.it.iastate.edu/accessing-clusters (VPN at vpn.iastate.edu)
- https://www.it.iastate.edu/remotework (Cisco Secure Client; links KB0011105)
